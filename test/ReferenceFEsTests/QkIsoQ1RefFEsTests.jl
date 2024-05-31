module QkIsoQ1RefFEsTests
  using Gridap
  using FillArrays
  using Test

  function generate_analytical_functions(order)
    u(x)   = x[1]^order
    f(x)   = -Δ(u)(x)
    u, f  
  end

  for D=1:3
    domain=Vector{Float64}(undef,D*2)
    partition=Vector{Int64}(undef,D)
    for d=1:D
      domain[2*d-1] = 0.0
      domain[2*d  ] = 1.0
      partition[d]  = 2
    end
    order  = 2
    I(x)   = x[1]^order
    u,f=generate_analytical_functions(order)

    modelH=CartesianDiscreteModel(Tuple(domain),Tuple(partition))
    modelh=Gridap.Adaptivity.refine(modelH,order) # For integration purposes

    reffeH = ReferenceFE(lagrangian,Float64,order)
    VH     = TestFESpace(modelH,reffeH;
                        conformity=:H1,
                        dirichlet_tags="boundary")
    UH     = TrialFESpace(VH,u)

    duH=get_trial_fe_basis(UH)

    reffeh=Gridap.QkIsoQ1(Float64,D,order)

    Vh=FESpace(modelH,reffeh;conformity=:H1,dirichlet_tags="boundary")
    Uh=TrialFESpace(Vh,u)

    ΩH   = Triangulation(modelH)
    Ωh   = Triangulation(modelh)
    dΩH  = Measure(ΩH,2*order+1)
    dΩh  = Measure(Ωh,2*order+1)
    dΩHh = Measure(ΩH,Ωh,2*order+1)

    ldΩHh(vh)=∫(1.0*vh)dΩHh
    b1ΩHh=assemble_vector(ldΩHh,Vh)

    ldΩH(vh)=∫(1.0*vh)dΩH
    b1ΩH=assemble_vector(ldΩH,Vh)

    println("norm(b1ΩHh-b1ΩH)=",norm(b1ΩHh-b1ΩH)) # High error with the "wrong" measure

    reffe = ReferenceFE(lagrangian,Float64,1)
    Vhl   = TestFESpace(modelh,reffe; conformity=:H1,dirichlet_tags="boundary")
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

    function generate_bilinear_form(dΩH)
      aH(u,v) = ∫(∇(u)⋅∇(v))*dΩH
    end 

    function generate_linear_form(f,dΩH)
      lH(v) = ∫(f*v)*dΩH
    end 

    function generate_weak_residual_form(Vh)
      function weak_res(uH, dΩ, f)
        vh_basis = Gridap.get_fe_basis(Vh)
        ∫(∇(uH)⋅∇(vh_basis))*dΩ - ∫(f*vh_basis)*dΩ
      end
    end

    for orderTrial in 1:2
      u,f = generate_analytical_functions(orderTrial)
      
      reffeH = ReferenceFE(lagrangian,Float64,orderTrial)
      VH = TestFESpace(modelH,reffeH; conformity=:H1, dirichlet_tags="boundary")
      UH = TrialFESpace(VH,u)
      ΩH = Triangulation(modelH)
      dΩH = Measure(ΩH,2*orderTrial+1)
      
      aH=generate_bilinear_form(dΩH)
      lH=generate_linear_form(f,dΩH)
      weak_res = generate_weak_residual_form(VH)
      
      # Test laplacian (Galerkin problem)
      # uH is the Galerkin solution
      op = AffineFEOperator(aH,lH,UH,VH)
      uH = solve(op)
      weak_res = generate_weak_residual_form(VH)
      rH = weak_res(uH,dΩH,f)
      ref_rnorm = norm(assemble_vector(rH,VH))
      @test ref_rnorm < 1.0e-12

      for orderLinearFERefine in 1:4
        reffeh = ReferenceFE(lagrangian,Float64,orderLinearFERefine)
        Vh = FESpace(modelH,reffeh;conformity=:H1,dirichlet_tags="boundary")
        order = max(orderTrial, orderLinearFERefine)
        Ωh = Triangulation(modelh)
        dΩHh = Measure(ΩH,Ωh,2*order+1)
        weak_res = generate_weak_residual_form(Vh)
        adΩHh(UH,vh)=∫(∇(UH)⋅∇(vh))dΩHh
        lΩHh(vh)    =∫(f*vh)dΩHh
        if orderTrial == orderLinearFERefine
          # Test laplacian (Petrov-Galerkin)
          op=AffineFEOperator(adΩHh,lΩHh,UH,Vh)
          uh=solve(op)
          eh=sum(∫((u-uh)*(u-uh))dΩHh)
          @test eh < 1.0e-12
          rh = weak_res(uh, dΩHh, f)
          rh_vec = assemble_vector(rh,Vh)
          @test norm(rh_vec) < 1.0e-10
        else 
          # Test Galerkin solution against LinearizedFESpace
          # 1. Using weak residual form 
          rh = weak_res(uH, dΩHh, f)
          rh_vec = assemble_vector(rh,Vh)
          @test norm(rh_vec) < 1.0e-06

          # 2. Using linear algebra data structures
          assem = SparseMatrixAssembler(UH,Vh)
          duH=get_trial_fe_basis(UH)
          dvh=get_fe_basis(Vh)
          data = Gridap.FESpaces.collect_cell_matrix_and_vector(UH,Vh,adΩHh(duH,dvh),lΩHh(dvh),zero(UH))
          A,b = assemble_matrix_and_vector(assem,data)
          x=get_free_dof_values(uH)
          r = A*x-b
          @test norm(r) < 1.0e-06
        end
      end
    end
  end   
end
 
