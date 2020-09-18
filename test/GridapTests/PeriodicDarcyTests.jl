module PeriodicDarcyTests

using Gridap
using Test
using LinearAlgebra

u(x) = VectorValue(x[1]*(x[1]-1)*(2x[2]-1.0),-x[2]*(x[2]-1.0)*(2x[1]-1.0))
p(x) = x[2]-0.5
f(x) = VectorValue(x[1]*(x[1]-1)*(2x[2]-1.0),-x[2]*(x[2]-1.0)*(2x[1]-1.0)+1.0)
g(x) = 0.0

domain = (0,1,0,1,0)
order = 1
partition = (4,4)
model = CartesianDiscreteModel(domain,partition; isperiodic=(true,false))

V = FESpace(
reffe=:RaviartThomas, order=order, valuetype=VectorValue{2,Float64},
conformity=:Hdiv, model=model,dirichlet_tags=collect(1:6))

Q = FESpace(
reffe=:QLagrangian, order=order, valuetype=Float64,
conformity=:L2, model=model, constraint=:zeromean)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

trian = Triangulation(model)
degree = 2*(order)
dΩ = LebesgueMeasure(trian,degree)
x = get_physical_coordinate(trian)

@form a((u,p),(v,q)) = ∫( u⋅v - p*(∇⋅v) + q*(∇⋅u) )*dΩ
@form l((v,q)) = ∫( v⋅f + q*(∇⋅u) )*dΩ
uh, ph = solve( a((U,P),(V,Q))==l((V,Q)) )

eu = u - uh
ep = p - ph

l2(v) = v⋅v
eu_l2 = sqrt(sum(∫(l2(eu))*dΩ))
ep_l2 = sqrt(sum(∫(l2(ep))*dΩ))

tol = 1.0e-9
@test eu_l2 < tol
@test ep_l2 < tol


end # module
