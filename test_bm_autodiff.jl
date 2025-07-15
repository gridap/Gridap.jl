using Gridap, Gridap.MultiField, Gridap.FESpaces
using BenchmarkTools, ProfileView, Profile

function benchmark(V,uh,res)
  ## Gradient
  println("Gradient:")
  _res_grad(uh) = res(uh,uh)
  _res_grafuh = _res_grad(uh)
  g = gradient(_res_grad,uh)
  G = assemble_vector(g,V)
  @btime Gridap.FESpaces._gradient($_res_grad,$uh,$_res_grafuh);
  @btime assemble_vector!($g,$G,$V);

  # ## Jacobian
  # println("Jacobian:")
  # dv = get_fe_basis(V)
  # _res_jac(uh) = res(uh,dv)
  # j = jacobian(_res_jac,uh)
  # J = assemble_matrix(j,V,V)
  # @btime jacobian($_res_jac,$uh);
  # @btime assemble_matrix!($j,$J,$V,$V);

  # ## Hessian
  # println("Hessian:")
  # h = hessian(_res_grad,uh);
  # H = assemble_matrix(h,V,V)
  # @btime hessian($_res_grad,$uh);
  # @btime assemble_matrix!($h,$H,$V,$V);
end

function benchmark_new(V,uh,res)
  ## Gradient
  println("Gradient:")
  _res_grad(uh) = res(uh,uh)
  _res_grafuh = _res_grad(uh)
  g = Gridap.MultiField._mf_gradient(f,uh,fuh)(_res_grad,uh,_res_grafuh)
  G = assemble_vector(g,V)
  @btime Gridap.MultiField._mf_gradient(f,uh,fuh)($_res_grad,$uh,$_res_grafuh);
  @btime assemble_vector!($g,$G,$V);

  # ## Jacobian
  # println("Jacobian:")
  # dv = get_fe_basis(V)
  # _res_jac(uh) = res(uh,dv)
  # j = jacobian(_res_jac,uh)
  # J = assemble_matrix(j,V,V)
  # @btime jacobian($_res_jac,$uh);
  # @btime assemble_matrix!($j,$J,$V,$V);

  # ## Hessian
  # println("Hessian:")
  # h = hessian(_res_grad,uh);
  # H = assemble_matrix(h,V,V)
  # @btime hessian($_res_grad,$uh);
  # @btime assemble_matrix!($h,$H,$V,$V);
end

function setup_mf(n,T,order)
  println("Running MF: n = $n, T = $T, order = $order")
  model = CartesianDiscreteModel((0,1,0,1),(n,n))
  Ω = Triangulation(model)
  dΩ = Measure(Ω,2*order)

  V = FESpace(Ω,ReferenceFE(lagrangian,T,order);dirichlet_tags="boundary")
  X = MultiFieldFESpace([V,V])
  uh = FEFunction(X,rand(num_free_dofs(X)))

  res((u1,u2),(v1,v2)) = ∫(∇(u1)⊙∇(v1) + (u1 ⋅ u2) * ∇(u2)⊙∇(v2))dΩ

  return X, uh, res
end

main_mf(n,T,order) = benchmark(setup_mf(n,T,order)...)
main_mf_new(n,T,order) = benchmark_new(setup_mf(n,T,order)...)

############################################################################################

main_mf(64,Float64,1);
main_mf_new(64,Float64,1);

main_mf(64,VectorValue{2,Float64},1);
main_mf_new(64,VectorValue{2,Float64},1);
