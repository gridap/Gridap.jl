module PoissonDGTests

using Test
using Gridap
import Gridap: ∇
#using LinearAlgebra

#domain = (0,1,0,1)
#partition = (4,4)
#model = CartesianDiscreteModel(domain,partition)
#const h = (domain[2]-domain[1]) / partition[1]

using Gridap.Geometry: DiscreteModelMock
model = DiscreteModelMock()
const h = 1

order = 2
const γ = 10

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)

degree = order
dΩ = LebesgueMeasure(Ω,degree)
dΓ = LebesgueMeasure(Γ,degree)
dΛ = LebesgueMeasure(Λ,degree)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)

u(x) = x[1]^2 + x[2]
∇u(x) = VectorValue( 2*x[1], one(x[2]) )
Δu(x) = 2
f(x) = - Δu(x)
∇(::typeof(u)) = ∇u

V = TestFESpace(model,ReferenceFE(:Lagrangian,Float64,order),conformity=:L2)
U = TrialFESpace(V,u)

a(u,v) =
  ∫( ∇(v)⋅∇(u) )*dΩ +
  ∫( (γ/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )*dΓ +
  ∫( (γ/h)*jump(v*n_Λ)⋅jump(u*n_Λ) - jump(v*n_Λ)⋅mean(∇(u)) -  mean(∇(v))⋅jump(u*n_Λ) )*dΛ

l(v) =
  ∫( v*f )*dΩ +
  ∫( (γ/h)*v*u - (n_Γ⋅∇(v))*u )*dΓ

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
