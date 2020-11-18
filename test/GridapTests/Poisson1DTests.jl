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
V = FESpace(model, ReferenceFE(:Lagrangian,Float64,order),conformity=:H1,dirichlet_tags=2)
U = TrialFESpace(V,u)

Ω = Triangulation(model) 
Γ = BoundaryTriangulation(model,1)
n_Γ = get_normal_vector(Γ)

degree = 2*order
dΩ = LebesgueMeasure(Ω,degree)
dΓ = LebesgueMeasure(Γ,degree)

a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
l(v) = ∫( v*f )*dΩ + ∫( v*(n_Γ⋅∇u) )*dΓ

op = AffineFEOperator(a,l,U,V)
uh = solve(op)

e = u - uh

l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
