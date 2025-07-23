using Gridap, Gridap.MultiField, Gridap.FESpaces
using BenchmarkTools

function benchmark(V,uh,res)
  ## Gradient
  println("Gradient:")
  _res_grad(uh) = res(uh,uh)
  g = gradient(_res_grad,uh;ad_type=:monolithic)
  G = assemble_vector(g,V)
  @btime gradient($_res_grad,$uh;ad_type=$(:monolithic));
  @btime assemble_vector!($g,$G,$V);

  ## Jacobian
  println("Jacobian:")
  dv = get_fe_basis(V)
  _res_jac(uh) = res(uh,dv)
  j = jacobian(_res_jac,uh;ad_type=:monolithic)
  J = assemble_matrix(j,V,V)
  @btime jacobian($_res_jac,$uh;ad_type=$(:monolithic));
  @btime assemble_matrix!($j,$J,$V,$V);

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
  g = gradient(_res_grad,uh;ad_type=:split)
  G = assemble_vector(g,V)
  @btime gradient($_res_grad,$uh;ad_type=$(:split));
  @btime assemble_vector!($g,$G,$V);

  ## Jacobian
  println("Jacobian:")
  dv = get_fe_basis(V)
  _res_jac(uh) = res(uh,dv)
  j = jacobian(_res_jac,uh;ad_type=:split)
  J = assemble_matrix(j,V,V)
  @btime jacobian($_res_jac,$uh;ad_type=$(:split));
  @btime assemble_matrix!($j,$J,$V,$V);

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

main_mf(64,VectorValue{2,Float64},2);
main_mf_new(64,VectorValue{2,Float64},2);