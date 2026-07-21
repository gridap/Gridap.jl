using Gridap

model = CartesianDiscreteModel((0,1,0,1), (3,3))
reffe = ReferenceFE(lagrangian, Float64, 1)
V = TestFESpace(model, reffe)
Ω = Skeleton(model)
dΩ = Measure(Ω, 2)
ener(u) = ∫( 1/2*jump(∇(u))⋅jump(∇(u)) )*dΩ
uh = interpolate(x->rand(),V)
hess = hessian(ener,uh)

du = get_trial_fe_basis(V)
dv = get_fe_basis(V)
hess_analytic = ∫( jump(∇(du))⋅jump(∇(dv)) )*dΩ

h = assemble_matrix(hess, V, V)
h_analytic = assemble_matrix(hess_analytic, V, V)

@show maximum(abs, h - h_analytic)

Gridap.Arrays.test_array(get_array(hess),get_array(hess_analytic),≈)