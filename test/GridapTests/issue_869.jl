module Issue869
using Gridap
using WriteVTK
using GridapEmbedded

d = mktempdir()

model = CartesianDiscreteModel(Point(0.0, 0.0), VectorValue(10.0, -1.0), (10, 5))

geo = disk(0.2, x0=0.5*VectorValue(10.0, -1.0))
cutgeo = cut(model, !(geo))
Ω = Interior(model)
Ω⁻ = Interior(cutgeo, PHYSICAL_IN)
Ω⁻act = Interior(cutgeo, ACTIVE_IN)
Γ = EmbeddedBoundary(cutgeo)

dΩ⁻act = Measure(Ω⁻act, 2)

dΩ⁻ = Measure(Ω⁻, 2)
dΓ = Measure(Γ, 2)

reffeᵤ = ReferenceFE(lagrangian, VectorValue{2, Float64}, 1)
reffeᵩ = ReferenceFE(lagrangian, Float64, 1)

Vstd = FESpace(Ω⁻act, reffeᵤ)
Ustd = TrialFESpace(Vstd)
Wstd = FESpace(Ω⁻act, reffeᵩ)
Φstd = TrialFESpace(Wstd)

X = MultiFieldFESpace([Ustd, Φstd])
Y = MultiFieldFESpace([Vstd, Wstd])

a((u, ϕ), (v, w)) = ∫(w*(u⋅VectorValue(1.0, 0.0)))dΓ                    

res((u, ϕ), (v, w)) = a((u, ϕ), (v, w))

jac(u, du, v) = jacobian(x->res(x,v),u)     
x = FEFunction(X, rand(num_free_dofs(X)))
dx = get_trial_fe_basis(X)
dy = get_fe_basis(Y)
jac(x,dx,dy)

op = FEOperator(res, jac, X, Y)

uh, ϕh = solve(op)
writevtk(Ω⁻act,joinpath(d,"tmp_model.vtu"),cellfields=["uh"=>uh, "phih"=>ϕh])
rm(d,recursive=true)
end