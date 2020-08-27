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
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model,1)
bquad = CellQuadrature(btrian,degree)
nb = get_normal_vector(btrian)

a(u,v) = ∇(v)⋅∇(u)
l(v) = v*f
t_Ω = AffineFETerm(a,l,trian,quad)

l_b(v) = v*(nb⋅∇u)
t_b = FESource(l_b,btrian,bquad)

op = AffineFEOperator(U,V,t_Ω,t_b)
uh = solve(op)

e = u - uh

l2(u) = inner(u,u)
sh1(u) = a(u,u)
h1(u) = sh1(u) + l2(u)

el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
ul2 = sqrt(sum( integrate(l2(uh),trian,quad) ))
uh1 = sqrt(sum( integrate(h1(uh),trian,quad) ))

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7


end # module
