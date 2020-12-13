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

V = FESpace(model,ReferenceFE(raviart_thomas,Float64,order),conformity=:Hdiv,
      dirichlet_tags=collect(1:6))


Q = FESpace(model,ReferenceFE(lagrangian,Float64,order);
      conformity=:L2, constraint=:zeromean)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

trian = Triangulation(model)
degree = 2*(order)
dΩ = Measure(trian,degree)
x = get_physical_coordinate(trian)

a((u, p),(v, q)) = ∫( u⋅v - (∇⋅v)*p + q*(∇⋅u) )*dΩ

b(( v, q)) = ∫( v⋅f + q*g)*dΩ


op = AffineFEOperator(a,b,X,Y)
xh = solve(op)
uh, ph = xh

eu = u - uh
ep = p - ph

l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
eu_l2 = l2(eu)
ep_l2 = l2(ep)

tol = 1.0e-9
@test eu_l2 < tol
@test ep_l2 < tol


end # module
