module LinearizedFESpacesTests
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
    modelH=CartesianDiscreteModel(Tuple(domain),Tuple(partition))

    order  = 4
    I(x)   = x[1]^order
    u,f=generate_analytical_functions(order)

    reffeH = ReferenceFE(lagrangian,Float64,order)
    VH     = TestFESpace(modelH,reffeH;
                        conformity=:H1,
                        dirichlet_tags="boundary")
    UH     = TrialFESpace(VH,u)


    Vh=Gridap.LinearizedFESpace(modelH,reffeH;
                                conformity=:H1,
                                dirichlet_tags="boundary")

    Uh=TrialFESpace(Vh,u)
    Ωh=get_triangulation(Vh)

    Ωh       = Triangulation(Vh.refined_model)
    dΩh      = Measure(Ωh,order+1)

    # Test L2-projection
    aL2dΩh(UH,vh)=∫(UH*vh)dΩh
    lL2Ωh(vh)    =∫(I*vh)dΩh
    op=AffineFEOperator(aL2dΩh,lL2Ωh,UH,Vh)
    Ih=solve(op)
    eh=sum(∫((I-Ih)*(I-Ih))dΩh)
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
        reffeh   = ReferenceFE(lagrangian,Float64,orderLinearFERefine)
        Vh       = Gridap.LinearizedFESpace(modelH,reffeh; conformity=:H1, dirichlet_tags="boundary")
        Ωh       = Triangulation(Vh.refined_model)
        dΩh      = Measure(Ωh,2*orderTrial+1)
        weak_res = generate_weak_residual_form(Vh)
        adΩHh(UH,vh)=∫(∇(UH)⋅∇(vh))dΩh
        lΩHh(vh)    =∫(f*vh)dΩh
        if orderTrial == orderLinearFERefine
          # Test laplacian (Petrov-Galerkin)
          op=AffineFEOperator(adΩHh,lΩHh,UH,Vh)
          uh=solve(op)
          eh=sum(∫((u-uh)*(u-uh))dΩh)
          @test eh < 1.0e-12
          rh = weak_res(uh, dΩh, f)
          rh_vec = assemble_vector(rh,Vh)
          @test norm(rh_vec) < 1.0e-12
        else
          # Test Galerkin solution against LinearizedFESpace
          # 1. Using weak residual form
          rh = weak_res(uH, dΩh, f)
          rh_vec = assemble_vector(rh,Vh)
          @test norm(rh_vec) < 1.0e-12

          # 2. Using linear algebra data structures
          assem = SparseMatrixAssembler(UH,Vh)
          duH=get_trial_fe_basis(UH)
          dvh=get_fe_basis(Vh)
          data = Gridap.FESpaces.collect_cell_matrix_and_vector(UH,Vh,adΩHh(duH,dvh),lΩHh(dvh),zero(UH))
          A,b = assemble_matrix_and_vector(assem,data)
          x=get_free_dof_values(uH)
          r = A*x-b
          @test norm(r) < 1.0e-12

          # 3. Boundary terms
          Γh = BoundaryTriangulation(Vh.refined_model, tags="boundary")
          dΓh = Measure(Γh, 2 * orderTrial + 1)
          n_Γh = get_normal_vector(Γh)
          @test length(assemble_vector(∫(n_Γh ⋅ ∇(dvh))dΓh, Vh)) > 0
        end
      end
    end
  end

  û(x) = sin(3.2 * x[1]^2) * cos(x[1]) + sin(4.6 * x[1]) * cos(5.2 * x[1])
  ŝ(x) = exp(x[1] / 2) + 2
  k̂(x) = 2 + sin(x[1])

  q(x) = k̂(x) * ∇(û)(x)
  f(x) = -(∇ ⋅ q)(x) + ŝ(x) * û(x)

  ncells = 10
  order = 5
  degree = 2 * order + 1
  reffe = ReferenceFE(lagrangian, Float64, order)
  model = CartesianDiscreteModel((0, 1), (ncells,))

  V0 = TestFESpace(model, reffe; dirichlet_tags="boundary")
  Ug = TrialFESpace(V0, û)
  Vl = Gridap.LinearizedFESpace(model, ReferenceFE(lagrangian, Float64, order); dirichlet_tags="boundary")
  Ω = Triangulation(Vl.refined_model)
  dΩ = Measure(Ω, degree)
  dv = Gridap.get_fe_basis(Vl)
  a(u, v) = ∫(∇(dv) ⋅ (k̂ * ∇(u)) + ŝ * u * dv)dΩ
  l(v) = ∫(f * dv)dΩ
  ũh = solve(AffineFEOperator(a, l, Ug, Vl))

  # jac = Gridap.jacobian(u -> a(u, dv) - l(dv), ũh)
  # A1=assemble_matrix(jac, Ug, Vl)
  # A2=assemble_matrix(a,Ug,Vl)
  # @assert norm(A1-A2) < 1.0e-12
end
