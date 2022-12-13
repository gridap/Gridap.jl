module LinearizedFESpacesTests
  using Gridap
  using FillArrays
  using Test

  for D=1:3
    domain=Vector{Float64}(undef,D*2)
    partition=Vector{Int64}(undef,D)
    for d=1:D
      domain[2*d-1] = 0.0
      domain[2*d  ] = 1.0
      partition[d]  = 2
    end
    modelH=CartesianDiscreteModel(Tuple(domain),Tuple(partition))

    order  = 2
    I(x)   = x[1]^order
    u(x)   = x[1]^order
    f(x)   = -Δ(u)(x)

    reffeH = ReferenceFE(lagrangian,Float64,order)
    VH     = TestFESpace(modelH,reffeH;
                        conformity=:H1,
                        dirichlet_tags="boundary")
    UH     = TrialFESpace(VH,u)

    duH=get_trial_fe_basis(UH)

    Vh=Gridap.LinearizedFESpace(modelH,reffeH;
                                conformity=:H1,
                                dirichlet_tags="boundary")

    Uh=TrialFESpace(Vh,u)

    ΩH   = Triangulation(modelH)
    Ωh   = Triangulation(Vh.refined_model)
    dΩH  = Measure(ΩH,2*order+1)
    dΩh  = Measure(Ωh,2*order+1)
    dΩHh = Measure(ΩH,Ωh,2*order+1)

    ldΩHh(vh)=∫(1.0*vh)dΩHh
    b1ΩHh=assemble_vector(ldΩHh,Vh)

    ldΩH(vh)=∫(1.0*vh)dΩH
    b1ΩH=assemble_vector(ldΩH,Vh)

    println("norm(b1ΩHh-b1ΩH)=",norm(b1ΩHh-b1ΩH)) # High error with the "wrong" measure

    reffe = ReferenceFE(lagrangian,Float64,1)
    Vhl   = TestFESpace(Vh.refined_model,reffe; conformity=:H1,dirichlet_tags="boundary")
    l(vh)=∫(1.0*vh)dΩh
    b2=assemble_vector(l,Vhl)

    @test norm(b1ΩHh-b2) < 1.0e-12 # Eps error in the "right" measure
    println("norm(b1ΩH-b2)=",norm(b1ΩH-b2)) # High error with the "wrong" measure

    # Test L2-projection
    aL2dΩHh(UH,vh)=∫(UH*vh)dΩHh
    lL2ΩHh(vh)    =∫(I*vh)dΩHh
    op=AffineFEOperator(aL2dΩHh,lL2ΩHh,UH,Vh)
    Ih=solve(op)
    eh=sum(∫((I-Ih)*(I-Ih))dΩHh)
    @test eh < 1.0e-12

    # Test laplacian
    ## VERY IMPORTANT: gradient of LinearizedFESpace basis MUST NOT BE computed
    #                  as gradient(vh). We need this specialized method.
    grad_vh=Gridap.get_gradient_fe_basis(Vh)
    adΩHh(UH,vh)=∫(∇(UH)⋅grad_vh)dΩHh
    lΩHh(vh)    =∫(f*vh)dΩHh
    op=AffineFEOperator(adΩHh,lΩHh,UH,Vh)
    uh=solve(op)
    eh=sum(∫((u-uh)*(u-uh))dΩHh)
    @test eh < 1.0e-12

    # Testing weak residual assembly when LinearizedFESpace is NOT
    # based on the refinement governed by TrialFESpace order, such
    # that the count and positions of DOFs in TrialFESpace and the
    # LinearizedFESpace do not remains the same

    aH(u,v) = ∫(∇(u)⋅∇(v))*dΩH
    lH(v) = ∫(f*v)*dΩH

    function weak_res(uH, Vh::Gridap.LinearizedFESpace, dΩ)
      grad_vh = Gridap.get_gradient_fe_basis(Vh)
      vh_basis = Gridap.get_fe_basis(Vh)
      ∫(∇(uH)⋅grad_vh)*dΩ - ∫(f*vh_basis)*dΩ
    end

    function weak_res(uH, VH, dΩ)
      vh_basis = Gridap.get_fe_basis(VH)
      ∫(∇(uH)⋅∇(vh_basis))*dΩ - ∫(f*vh_basis)*dΩ
    end

    for orderTrial in 1:2
      reffeH = ReferenceFE(lagrangian,Float64,orderTrial)
      VH = TestFESpace(modelH,reffeH; conformity=:H1, dirichlet_tags="boundary")
      UH = TrialFESpace(VH,u)
      ΩH = Triangulation(modelH)
      dΩH = Measure(ΩH,2*orderTrial+1)

      op = AffineFEOperator(aH,lH,UH,VH)
      uH = solve(op)
      rH = weak_res(uH, VH, dΩH)
      ref_rnorm = norm(assemble_vector(rH,VH))

      for orderLinearFERefine in 1:4
        reffeh = ReferenceFE(lagrangian,Float64,orderLinearFERefine)
        Vh = Gridap.LinearizedFESpace(modelH,reffeh; conformity=:H1, dirichlet_tags="boundary")
        order = max(orderTrial, orderLinearFERefine)
        Ωh = Triangulation(Vh.refined_model)
        dΩh = Measure(Ωh,2*orderLinearFERefine+1)
        dΩHh = Measure(ΩH,Ωh,2*order+1)
        rh = weak_res(uH, Vh, dΩHh)
        rh_vec = assemble_vector(rh,Vh)
        @test length(rh_vec) == num_free_dofs(Vh)
      end
    end

  end
end
