module Poisson1DTests

using Gridap
using Test

u(x) = x[1]^2 + 1
∇u(x) = ∇(u)(x)
f(x) = -Δ(u)(x)

domain = (0,3)
cells = (10,)
model = CartesianDiscreteModel(domain,cells)

order = 2
V = FESpace(
  model=model,
  reffe=:Lagrangian,
  order=order,
  valuetype=Float64,
  conformity=:H1,
  dirichlet_tags=2)

U = TrialFESpace(V,u)

degree = 2*order
trian = Triangulation(model) 
dΩ = LebesgueMeasure(trian,degree)

btrian = BoundaryTriangulation(model,1)
dΓ = LebesgueMeasure(btrian,degree)
nb = get_normal_vector(btrian)

@form a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
@form l(v) = ∫( v*f )*dΩ + ∫( v*(nb⋅∇u) )*dΓ

uh = solve( a(U,V)==l(V) )

e = u - uh

l2(u) = u⊙u
sh1(u) = ∇(u)⊙∇(u)
h1(u) = sh1(u) + l2(u)

el2 = sqrt(sum( ∫( l2(e) )*dΩ ))
eh1 = sqrt(sum( ∫( h1(e) )*dΩ ))
ul2 = sqrt(sum( ∫( l2(uh) )*dΩ ))
uh1 = sqrt(sum( ∫( h1(uh) )*dΩ ))

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7


end # module
