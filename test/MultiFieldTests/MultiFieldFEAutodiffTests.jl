module MultiFieldFEAutodiffTests

using Test

using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

Ω = Triangulation(model)
dΩ = Measure(Ω,2)

V1 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
V2 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
Y = MultiFieldFESpace([V1,V2])

r((uh,ph),(v,q)) = ∫( v*uh + uh*q + ph*q )dΩ
j(xh,dx,dy) = jacobian(xh->r(xh,dy),xh)

e((uh,ph)) = ∫( uh*uh + uh*ph + ph*ph )dΩ
g(xh,dy) = gradient(xh->e(xh),xh)
#h(xh,dx,dy) = hessian(xh->e(xh),xh)

dx = get_trial_fe_basis(Y)
dy = get_fe_basis(Y)
xh = FEFunction(Y,rand(num_free_dofs(Y)))

#display(r(xh,dy)[Ω][end])
#display(j(xh,dx,dy)[Ω][end])
#display(e(xh)[Ω][end])
#display(g(xh,dy)[Ω][end])
##display(h(xh,dx,dy)[Ω][end])

@test j(xh,dx,dy)[Ω][end][1,1] != nothing
@test j(xh,dx,dy)[Ω][end][2,1] != nothing
@test j(xh,dx,dy)[Ω][end][1,2] == nothing
@test j(xh,dx,dy)[Ω][end][2,2] != nothing
@test g(xh,dy)[Ω][end][1] != nothing
@test g(xh,dy)[Ω][end][2] != nothing

V1 = FESpace(model,ReferenceFE(lagrangian,Float64,2))
V2 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
Y = MultiFieldFESpace([V1,V2])

dx = get_trial_fe_basis(Y)
dy = get_fe_basis(Y)
xh = FEFunction(Y,rand(num_free_dofs(Y)))

@test_broken j(xh,dx,dy)[Ω][end][1,1] != nothing
@test_broken j(xh,dx,dy)[Ω][end][2,1] != nothing
@test_broken j(xh,dx,dy)[Ω][end][1,2] == nothing
@test_broken j(xh,dx,dy)[Ω][end][2,2] != nothing
@test_broken g(xh,dy)[Ω][end][1] != nothing
@test_broken g(xh,dy)[Ω][end][2] != nothing

end # module
