module StokesDGTests

using Test
using Gridap
import Gridap: ∇

u(x) = VectorValue(x[1]*x[1], x[2])
∇u(x) = TensorValue(2*x[1],0.0,0.0,1.0)
Δu(x) = VectorValue(2.0,0.0)

p(x) = x[1] - x[2]
∇p(x) = VectorValue(1.0,-1.0)

f(x) = - Δu(x) + ∇p(x)
g(x) = tr(∇u(x))

∇(::typeof(u)) = ∇u
∇(::typeof(p)) = ∇p

L = 1.0
domain = (0.0, L, 0.0, L)
ncellx = 6
partition = (ncellx,ncellx)
model = CartesianDiscreteModel(domain,partition)
model = simplexify(model)

order = 2
const h = L / ncellx
const γ = order*(order+1)
const γ0 = 1.0/10.0

reffe_u = ReferenceFE(:Lagrangian,VectorValue{2,Float64},order)
reffe_p = ReferenceFE(:Lagrangian,Float64,order)

U = FESpace(model,reffe_u,conformity=:L2)
P = FESpace(model,reffe_p,conformity=:L2,constraint=:zeromean)
X = MultiFieldFESpace([U,P])

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)

degree = 2*order
dΩ = LebesgueMeasure(Ω,degree)
dΓ = LebesgueMeasure(Γ,degree)
dΛ = LebesgueMeasure(Λ,degree)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)

a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) - ∇(q)⋅u + v⋅∇(p) )*dΩ +
  ∫( (γ/h)*v⋅u - v⋅(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))⋅u + 2*(q*n_Γ)⋅u )*dΓ +
  ∫(  
    (γ/h)*jump(v⊗n_Λ)⊙jump(u⊗n_Λ) -
      jump(v⊗n_Λ)⊙mean(∇(u)) -
      mean(∇(v))⊙jump(u⊗n_Λ)  +
      (γ0*h)*jump(q*n_Λ)⋅jump(p*n_Λ) +
      jump(q*n_Λ)⋅mean(u) -
      mean(v)⋅jump(p*n_Λ)
   )*dΛ

l((v,q)) =
  ∫( v⋅f + q*g )*dΩ +
  ∫( (γ/h)*v⋅u - (n_Γ⋅∇(v))⋅u + (q*n_Γ)⋅u )*dΓ

op = AffineFEOperator(a,l,X,X)

uh, ph = solve(op)

eu = u - uh
ep = p - ph

l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

eu_l2 = l2(eu)
eu_h1 = h1(eu)
ep_l2 = l2(ep)

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
