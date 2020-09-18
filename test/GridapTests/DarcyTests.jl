module DarcyTests

using Test
using Gridap
import Gridap: ∇, divergence
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

V = FESpace(
  reffe=:RaviartThomas, order=order, valuetype=VectorValue{2,Float64},
  conformity=:Hdiv, model=model, dirichlet_tags=[5,6])

Q = FESpace(
  reffe=:QLagrangian, order=order, valuetype=Float64,
  conformity=:L2, model=model)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

trian = Triangulation(model)
degree = 2
dΩ = LebesgueMeasure(trian,degree)

neumanntags = [7,8]
btrian = BoundaryTriangulation(model,neumanntags)
degree = 2*(order+1)
dΓ = LebesgueMeasure(btrian,degree)
nb = get_normal_vector(btrian)

@form a((u,p),(v,q)) = 
  ∫( u⋅v - p*(∇⋅v) + q*(∇⋅u) )*dΩ

@form l((v,q)) = 
  ∫( v⋅f + q*(∇⋅u) )*dΩ +
  ∫( -(v⋅nb)*p )*dΓ

uh, ph = solve( a((U,P),(V,Q))==l((V,Q)) )

eu = u - uh
ep = p - ph

l2(v) = v⋅v
h1(v) = v⋅v + ∇(v)⊙∇(v)

eu_l2 = sqrt(sum(∫(l2(eu))*dΩ))
eu_h1 = sqrt(sum(∫(h1(eu))*dΩ))
ep_l2 = sqrt(sum(∫(l2(ep))*dΩ))

tol = 1.0e-8
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
