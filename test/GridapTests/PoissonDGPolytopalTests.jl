module PoissonDGPolytopalTests

using Test
using Gridap
using LinearAlgebra
using Gridap.FESpaces

D = 3
domain = Tuple(repeat([0, 1], D))
partition = Tuple(fill(4, D))
base_model = simplexify(CartesianDiscreteModel(domain,partition); positive=true)

if D == 2
  model = Gridap.Geometry.voronoi(base_model)
else
  model = Gridap.Geometry.PolytopalDiscreteModel(base_model)
end

order = 2
γ = 10

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)

degree = order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΛ = Measure(Λ,degree)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)
h_Γ = CellField(1 ./ get_array(∫(1)*dΓ) .^(D-1), Γ)
h_Λ = CellField(1 ./ get_array(∫(1)*dΛ) .^(D-1), Λ) 

u(x) = x[1]^2 + x[2]
f(x) = - Δ(u)(x)

V = FESpaces.PolytopalFESpace(model,Float64,order)

a(u,v) =
  ∫( ∇(v)⋅∇(u) )*dΩ +
  ∫( (γ/h_Γ)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )*dΓ +
  ∫( (γ/h_Λ)*jump(v*n_Λ)⋅jump(u*n_Λ) - jump(v*n_Λ)⋅mean(∇(u)) -  mean(∇(v))⋅jump(u*n_Λ) )*dΛ

l(v) =
  ∫( v*f )*dΩ +
  ∫( (γ/h_Γ)*v*u - (n_Γ⋅∇(v))*u )*dΓ

op = AffineFEOperator(a,l,V,V)
uh = solve(op)

eh = u - uh
el2 = sqrt(sum( ∫( eh⊙eh )*dΩ ))
eh1 = sqrt(sum( ∫( eh⊙eh + ∇(eh)⊙∇(eh) )*dΩ ))

@test el2 < 1.e-8
@test eh1 < 1.e-7

end # module
