module MultiFieldFEAutodiffTests

using Test

using Gridap
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

using ForwardDiff
using SparseArrays

function run_tests(ad_type)
  domain = (0,1,0,1)
  partition = (2,2)
  model = CartesianDiscreteModel(domain,partition)

  Ω = Triangulation(model)
  dΩ = Measure(Ω,2)

  V1 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
  V2 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
  Y = MultiFieldFESpace([V1,V2])

  r((uh,ph),(v,q)) = ∫( v*uh + uh*q + ph*q )dΩ
  j(xh,dx,dy) = jacobian(xh->r(xh,dy),xh;ad_type)

  e((uh,ph)) = ∫( uh*uh + uh*ph + ph*ph )dΩ
  g(xh,dy) = gradient(xh->e(xh),xh;ad_type)
  h(xh,dx,dy) = hessian(xh->e(xh),xh;ad_type)

  dx = get_trial_fe_basis(Y)
  dy = get_fe_basis(Y)
  xh = FEFunction(Y,rand(num_free_dofs(Y)))

  #display(r(xh,dy)[Ω][end])
  #display(j(xh,dx,dy)[Ω][end])
  #display(e(xh)[Ω][end])
  #display(g(xh,dy)[Ω][end])
  ##display(h(xh,dx,dy)[Ω][end])

  @test !isnothing(j(xh,dx,dy)[Ω][end][1,1])
  @test !isnothing(j(xh,dx,dy)[Ω][end][2,1])
  @test !isnothing(j(xh,dx,dy)[Ω][end][1,2])
  @test !isnothing(j(xh,dx,dy)[Ω][end][2,2])
  @test !isnothing(g(xh,dy)[Ω][end][1])
  @test !isnothing(g(xh,dy)[Ω][end][2])
  if ad_type == :monolithic
    @test !isnothing(h(xh,dx,dy)[Ω][end][1,1])
    @test !isnothing(h(xh,dx,dy)[Ω][end][2,1])
    @test !isnothing(h(xh,dx,dy)[Ω][end][1,2])
    @test !isnothing(h(xh,dx,dy)[Ω][end][2,2])
  end

  V1 = FESpace(model,ReferenceFE(lagrangian,Float64,2))
  V2 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
  Y = MultiFieldFESpace([V1,V2])

  dx = get_trial_fe_basis(Y)
  dy = get_fe_basis(Y)
  xh = FEFunction(Y,rand(num_free_dofs(Y)))

  @test !isnothing(j(xh,dx,dy)[Ω][end][1,1])
  @test !isnothing(j(xh,dx,dy)[Ω][end][2,1])
  @test !isnothing(j(xh,dx,dy)[Ω][end][1,2])
  @test !isnothing(j(xh,dx,dy)[Ω][end][2,2])
  @test !isnothing(g(xh,dy)[Ω][end][1])
  @test !isnothing(g(xh,dy)[Ω][end][2])
  if ad_type == :monolithic
    @test !isnothing(h(xh,dx,dy)[Ω][end][1,1])
    @test !isnothing(h(xh,dx,dy)[Ω][end][2,1])
    @test !isnothing(h(xh,dx,dy)[Ω][end][1,2])
    @test !isnothing(h(xh,dx,dy)[Ω][end][2,2])
  end

  eu(uh) = ∫( uh*uh )dΩ
  ep(ph) = ∫( ph*ph )dΩ
  ez((uh,ph)) = ∫( 0*uh*ph )dΩ
  eup((uh,ph)) = ∫( uh*uh + ph*ph )dΩ
  geu(xh,dy) = gradient(xh->eu(xh),xh)
  gep(xh,dy) = gradient(xh->ep(xh),xh)
  geup(xh,dy) = gradient(xh->eup(xh),xh;ad_type)
  heu(xh,dx,dy) = hessian(xh->eu(xh),xh)
  hep(xh,dx,dy) = hessian(xh->ep(xh),xh)
  hez(xh,dx,dy) = hessian(xh->ez(xh),xh)
  heup(xh,dx,dy) = hessian(xh->eup(xh),xh;ad_type)

  uh,ph=xh
  geu_uh=geu(uh,dy)
  gep_ph=gep(ph,dy)
  geup_uh_ph=geup(xh,dy)

  a=geu_uh[Ω]
  b=gep_ph[Ω]
  c=geup_uh_ph[Ω]

  @test all(a .== map(x->x.array[1],c))
  @test all(b .== map(x->x.array[2],c))

  if ad_type==:monolithic
    heu_uh=heu(uh,dy,dy)
    hep_ph=hep(ph,dy,dy)
    hez_uh=hez(xh,dy,dy)
    heup_uh_ph=heup(xh,dy,dy)

    a=heu_uh[Ω]
    b=hep_ph[Ω]
    c=hez_uh[Ω]
    d=heup_uh_ph[Ω]

    @test all(a .== map(x->x.array[1,1],d))
    @test all(b .== map(x->x.array[2,2],d))
    @test all(map(x->x.array[2,1],c) .== map(x->x.array[2,1],d))
    @test all(map(x->x.array[1,2],c) .== map(x->x.array[1,2],d))
  end

  ru(uh,v) = ∫( v*uh*uh )dΩ
  rp(ph,q) = ∫( q*ph*ph )dΩ
  rup((uh,ph),(v,q)) = ∫( v*uh*uh + q*ph*ph )dΩ
  ju(xh,dx,dy) = jacobian(xh->ru(xh,dy),xh)
  jp(xh,dx,dy) = jacobian(xh->rp(xh,dy),xh)
  jup(xh,dx,dy) = jacobian(xh->rup(xh,dy),xh;ad_type)

  uh=FEFunction(V1,get_free_dof_values(xh[1]))
  ph=FEFunction(V2,get_free_dof_values(xh[2]))
  du = get_trial_fe_basis(V1)
  dv = get_fe_basis(V1)
  dp = get_trial_fe_basis(V2)
  dq = get_fe_basis(V2)


  ju_uh=ju(uh,du,dv)
  jp_ph=jp(ph,dp,dq)
  jup_uh_ph=jup(xh,dx,dy)

  a=ju_uh[Ω]
  b=jp_ph[Ω]
  c=jup_uh_ph[Ω]

  @test all(a .== map(x->x.array[1,1],c))
  @test all(b .== map(x->x.array[2,2],c))

  # AD tests for Multifiled functionals involving SkeletonTriangulation
  # compared with direct ForwardDiff results

  reffeV = ReferenceFE(lagrangian,Float64,2)
  reffeQ = ReferenceFE(lagrangian,Float64,1)
  V = TestFESpace(model,reffeV,conformity=:L2)
  Q = TestFESpace(model,reffeQ,conformity=:L2)

  U = TrialFESpace(V)
  P = TrialFESpace(Q)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  xh = FEFunction(X,rand(num_free_dofs(X)))

  cell_u = get_cell_dof_values(xh)
  xh_bis = CellField(X,cell_u)
  cf = SkeletonCellFieldPair(xh,xh_bis)

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ,2)
  n_Λ = get_normal_vector(Λ)

  g_Λ((uh,ph)) = ∫( mean(uh) + mean(ph) + mean(uh)*mean(ph) )dΛ
  f_Λ((uh,ph)) = ∫( mean(uh*uh) + mean(uh*ph) + mean(ph*ph) )dΛ
  a_Λ((uh,ph)) = ∫( - jump(uh*n_Λ)⊙mean(∇(ph))
                    - mean(∇(uh))⊙jump(ph*n_Λ)
                    + jump(uh*n_Λ)⊙jump(ph*n_Λ) )dΛ

  function f_uh_free_dofs(f,xh,θ)
    dir_u = similar(xh[1].dirichlet_values,eltype(θ))
    dir_p = similar(xh[2].dirichlet_values,eltype(θ))
    X = xh.fe_space
    θ_uh = restrict_to_field(X,θ,1)
    θ_ph = restrict_to_field(X,θ,2)
    uh = FEFunction(X[1],θ_uh,dir_u)
    ph = FEFunction(X[2],θ_ph,dir_p)
    xh = MultiFieldFEFunction(θ,xh.fe_space,[uh,ph])
    sum(f(xh))
  end
  # can also do the above by constructing a MultiFieldCellField

  g_Λ_(θ) = f_uh_free_dofs(g_Λ,xh,θ)
  f_Λ_(θ) = f_uh_free_dofs(f_Λ,xh,θ)
  a_Λ_(θ) = f_uh_free_dofs(a_Λ,xh,θ)

  θ = get_free_dof_values(xh)

  # check if the evaluations are working
  @test sum(f_Λ(xh)) == f_Λ_(θ)
  @test sum(a_Λ(xh)) == a_Λ_(θ)
  @test sum(g_Λ(xh)) == g_Λ_(θ)

  f_Ω((uh,ph)) = ∫(uh*uh + ph*uh + uh*ph)*dΩ
  f_Ω_(θ) = f_uh_free_dofs(f_Ω,xh,θ)
  gridapgradf_Ω = assemble_vector(gradient(f_Ω,xh;ad_type),X)
  fdgradf_Ω = ForwardDiff.gradient(f_Ω_,θ)
  test_array(gridapgradf_Ω,fdgradf_Ω,≈)

  gridapgradf = assemble_vector(gradient(f_Λ,xh;ad_type),X)
  fdgradf = ForwardDiff.gradient(f_Λ_,θ)
  test_array(gridapgradf,fdgradf,≈)

  gridapgradg = assemble_vector(gradient(g_Λ,xh;ad_type),X)
  fdgradg = ForwardDiff.gradient(g_Λ_,θ)
  test_array(gridapgradg,fdgradg,≈)

  gridapgrada = assemble_vector(gradient(a_Λ,xh;ad_type),X)
  fdgrada = ForwardDiff.gradient(a_Λ_,θ)
  test_array(gridapgrada,fdgrada,≈)

  function jac_uh_free_dofs(j,xh,θ)
    dir_u = similar(xh[1].dirichlet_values,eltype(θ))
    dir_p = similar(xh[2].dirichlet_values,eltype(θ))
    X = xh.fe_space
    θ_uh = restrict_to_field(X,θ,1)
    θ_ph = restrict_to_field(X,θ,2)
    uh = FEFunction(X[1],θ_uh,dir_u)
    ph = FEFunction(X[2],θ_ph,dir_p)
    xh = MultiFieldFEFunction(θ,xh.fe_space,[uh,ph])

    T = eltype(θ)
    matrix_type = SparseMatrixCSC{T,Int}
    vector_type = Vector{T}
    assem = SparseMatrixAssembler(matrix_type,vector_type,X,X)
    assemble_vector(v->j(xh,v), assem, X)
  end

  J_Λ((uh,ph),(duh,dph)) = ∫( mean(duh)*mean(uh)*mean(uh) + mean(dph)*mean(ph)*mean(uh) + mean(duh)*mean(ph) )dΛ
  J_Λ_(θ) = jac_uh_free_dofs(J_Λ,xh,θ)

  gridapjacf = assemble_matrix((du,dv) -> jacobian(u -> J_Λ(u,dv),xh;ad_type),X,X)
  fdjacf = ForwardDiff.jacobian(J_Λ_,θ)
  test_array(gridapjacf,fdjacf,≈)

end

function run_multi_trian_tests()
  domain = (0,1,0,1)
  partition = (2,2)
  model = CartesianDiscreteModel(domain,partition)
  Ω = Triangulation(model)
  Γ = BoundaryTriangulation(model,tags=["tag_5"])
  dΩ = Measure(Ω,2)
  dΓ = Measure(Γ,2)
  V1 = FESpace(Γ,ReferenceFE(lagrangian,Float64,1))
  V2 = FESpace(model,ReferenceFE(lagrangian,VectorValue{2,Float64},1))
  V3 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
  X = MultiFieldFESpace([V1,V2,V3])
  f(xh) = ∫(xh[1]+xh[2]⋅xh[2]+xh[1]*xh[3])dΓ
  uh = zero(X)
  du = gradient(f,uh;ad_type=:split)
  du1 = gradient(x->f((x,uh[2],uh[3])),uh[1])
  du2 = gradient(x->f((uh[1],x,uh[3])),uh[2])
  du3 = gradient(x->f((uh[1],uh[2],x)),uh[3])

  @test lazy_map(Base.Fix2(getindex,1),du[Γ]) == du1[Γ]
  @test lazy_map(Base.Fix2(getindex,2),du[Γ]) == du2[Γ]
  @test lazy_map(Base.Fix2(getindex,3),du[Γ]) == du3[Γ]

  du1_vec = assemble_vector(du1,V1)
  du2_vec = assemble_vector(du2,V2)
  du3_vec = assemble_vector(du3,V3)
  du_vec = assemble_vector(du,X)

  @test du_vec == [du1_vec;du2_vec;du3_vec]

  f2(xh,yh) = ∫(xh[1]⋅yh[1]+xh[2]⋅yh[2]+xh[1]⋅xh[2]⋅yh[2]+xh[1]*xh[3]*yh[3])dΓ
  dv = get_fe_basis(X)
  _j = jacobian(uh->f2(uh,dv),uh;ad_type=:split)
  J = assemble_matrix(_j,X,X)

  f2_jac(xh,dxh,yh) = ∫(dxh[1]⋅yh[1]+dxh[2]⋅yh[2]+dxh[1]⋅xh[2]⋅yh[2]+xh[1]⋅dxh[2]⋅yh[2]+dxh[1]*xh[3]*yh[3]+xh[1]*dxh[3]*yh[3])dΓ
  op = FEOperator(f2,f2_jac,X,X)
  J_fwd = jacobian(op,uh)

  @test J_fwd == J

  V1 = FESpace(Triangulation(model,1:2),ReferenceFE(lagrangian,Float64,1),conformity=:L2)
  V2 = FESpace(model,ReferenceFE(lagrangian,VectorValue{2,Float64},1),conformity=:L2)
  V3 = FESpace(model,ReferenceFE(lagrangian,Float64,1),conformity=:L2)
  X = MultiFieldFESpace([V1,V2,V3])
  uh = zero(X)

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ,2)
  _f(xh) = ∫(mean(xh[1])+mean(xh[2])⋅mean(xh[2])+mean(xh[1])*mean(xh[3]))dΛ
  uh = zero(X)
  du = gradient(_f,uh;ad_type=:split)
  du1 = gradient(x->_f((x,uh[2],uh[3])),uh[1])
  du2 = gradient(x->_f((uh[1],x,uh[3])),uh[2])
  du3 = gradient(x->_f((uh[1],uh[2],x)),uh[3])

  @test lazy_map(Base.Fix2(getindex,1),du[Λ]) == du1[Λ]
  @test lazy_map(Base.Fix2(getindex,2),du[Λ]) == du2[Λ]
  @test lazy_map(Base.Fix2(getindex,3),du[Λ]) == du3[Λ]

  du1_vec = assemble_vector(du1,V1)
  du2_vec = assemble_vector(du2,V2)
  du3_vec = assemble_vector(du3,V3)
  du_vec = assemble_vector(du,X)

  @test du_vec == [du1_vec;du2_vec;du3_vec]

  _f2(xh,yh) = ∫(mean(xh[1])⋅mean(yh[1])+mean(xh[2])⋅mean(yh[2])+mean(xh[1])⋅mean(xh[2])⋅mean(yh[2])+mean(xh[1])*mean(xh[3])*mean(yh[3]))dΛ
  dv = get_fe_basis(X)
  _j = jacobian(uh->_f2(uh,dv),uh;ad_type=:split)
  J = assemble_matrix(_j,X,X)

  _f2_jac(xh,dxh,yh) = ∫(mean(dxh[1])⋅mean(yh[1])+mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])⋅mean(xh[2])⋅mean(yh[2]) +
    mean(xh[1])⋅mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])*mean(xh[3])*mean(yh[3])+mean(xh[1])*mean(dxh[3])*mean(yh[3]))dΛ
  op = FEOperator(_f2,_f2_jac,X,X)
  J_fwd = jacobian(op,uh)

  @test J_fwd == J
end

function run_edge_case_test()
  for mod_type in (:skel_trian, :trian, :boundary)
    for style in (ConsecutiveMultiFieldStyle(),BlockMultiFieldStyle())
      model = CartesianDiscreteModel((0,1,0,1),(3,3))
      Ω1 = Triangulation(model,[1,2,3])
      if mod_type == :skel_trian
        Λ = SkeletonTriangulation(model)
      elseif mod_type == :trian
        Λ = Triangulation(model)
      elseif mod_type == :boundary
        Λ = BoundaryTriangulation(model)
      else
        error()
      end

      V1 = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1),conformity=:L2)
      V2 = FESpace(Ω1,ReferenceFE(lagrangian,VectorValue{2,Float64},1),conformity=:L2)
      V3 = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1),conformity=:L2)
      V = MultiFieldFESpace([V1,V2,V3];style)
      uh = zero(V);
      fill!(get_free_dof_values(uh),1.0)
      dΛ = Measure(Λ,2)

      # Gradient
      f(xh) = ∫(mean(xh[1])+mean(xh[2])⋅mean(xh[2])+mean(xh[1])*mean(xh[3]))dΛ
      dv = get_fe_basis(V);
      j = gradient(f,uh;ad_type=:split)
      J = assemble_vector(j,V)

      j_mono = gradient(f,uh;ad_type=:monolithic)
      J_mono = assemble_vector(j_mono,V)

      df2(v,xh) = ∫(mean(v[1])+2*mean(v[2])⋅mean(xh[2])+mean(v[1])*mean(xh[3])+mean(xh[1])*mean(v[3]))dΛ
      J_analytic = assemble_vector(dv->df2(dv,uh),V)
      @test J ≈ J_analytic
      @test J_mono ≈ J_analytic

      # Jacobian
      f2(xh,yh) = ∫(mean(xh[1])⋅mean(yh[1])+mean(xh[2])⋅mean(yh[2])+mean(xh[1])⋅mean(xh[2])⋅mean(yh[2])+mean(xh[1])*mean(xh[3])*mean(yh[3]))dΛ
      dv = get_fe_basis(V);
      j = jacobian(uh->f2(uh,dv),uh;ad_type=:split)
      J = assemble_matrix(j,V,V)

      j_mono = jacobian(uh->f2(uh,dv),uh;ad_type=:monolithic)
      J_mono = assemble_matrix(j_mono,V,V)

      df2(xh,dxh,yh) = ∫(mean(dxh[1])⋅mean(yh[1])+mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])⋅mean(xh[2])⋅mean(yh[2]) +
        mean(xh[1])⋅mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])*mean(xh[3])*mean(yh[3])+mean(xh[1])*mean(dxh[3])*mean(yh[3]))dΛ
      op = FEOperator(f2,df2,V,V)
      J_analytic = jacobian(op,uh)

      @test J ≈ J_analytic
    end
  end
end

run_tests(:split)
run_tests(:monolithic)
run_multi_trian_tests()
run_edge_case_test()

end # module
