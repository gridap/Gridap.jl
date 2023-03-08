module Issue869
using Gridap

model = CartesianDiscreteModel(Point(0.0, 0.0), VectorValue(10.0, -1.0), (10, 5))

Ω = Interior(model)
mask = map(isodd, 1:num_cells(model))
Ω⁻act = Interior(model, mask)

dΩ = Measure(Ω, 2)

reffeᵤ = ReferenceFE(lagrangian, VectorValue{2, Float64}, 1)
reffeᵩ = ReferenceFE(lagrangian, Float64, 1)

Vstd = FESpace(Ω⁻act, reffeᵤ)
Ustd = TrialFESpace(Vstd)
Wstd = FESpace(Ω⁻act, reffeᵩ)
Φstd = TrialFESpace(Wstd)

X = MultiFieldFESpace([Ustd, Φstd])
Y = MultiFieldFESpace([Vstd, Wstd])

a((u, ϕ), (v, w)) = ∫(w*(u⋅VectorValue(1.0, 0.0)))dΩ                    

res((u, ϕ), (v, w)) = a((u, ϕ), (v, w))

jac(u, du, v) = jacobian(x->res(x, v), u)     
x = FEFunction(X, rand(num_free_dofs(X)))
dx = get_trial_fe_basis(X)
dy = get_fe_basis(Y)
jac(x, dx, dy)

end