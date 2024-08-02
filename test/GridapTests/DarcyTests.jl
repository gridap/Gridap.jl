module DarcyTests

using Test
using Gridap
import Gridap: ∇, divergence, DIV
using LinearAlgebra

u(x) = VectorValue(2*x[1],x[1]+x[2])

divergence(::typeof(u)) = (x) -> 3

p(x) = x[1]-x[2]

∇p(x) = VectorValue(1,-1)

∇(::typeof(p)) = ∇p

f(x) = u(x) + ∇p(x)

domain = (0,1,0,1)
partition = (4,4)
order = 1
model = CartesianDiscreteModel(domain,partition)

V = FESpace(model,ReferenceFE(raviart_thomas,Float64,order),conformity=:Hdiv,
  dirichlet_tags=[5,6])

Q = FESpace(model,ReferenceFE(lagrangian,Float64,order); conformity=:L2)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

trian = Triangulation(model)
degree = 2
dΩ = Measure(trian,degree)
dω = Measure(trian,degree,integration_domain_style=ReferenceDomain())

neumanntags = [7,8]
btrian = BoundaryTriangulation(model,tags=neumanntags)
degree = 2*(order+1)
dΓ = Measure(btrian,degree)
nb = get_normal_vector(btrian)

a((u, p),(v, q)) = ∫( u⋅v )*dΩ + ∫(q*DIV(u)-DIV(v)*p)*dω

b(( v, q)) = ∫( v⋅f + q*(∇⋅u))*dΩ - ∫((v⋅nb)*p )*dΓ

op = AffineFEOperator(a,b,X,Y)
xh = solve(op)
uh, ph = xh

eu = u - uh
ep = p - ph

l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
h1(v) = sqrt(sum(∫(v*v + ∇(v)⋅∇(v))*dΩ))

eu_l2 = l2(eu)
ep_l2 = l2(ep)
ep_h1 = h1(ep)

tol = 1.0e-9
@test eu_l2 < tol
@test ep_l2 < tol
@test ep_h1 < tol

end # module
