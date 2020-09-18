module StokesDGTests

using Test
using Gridap
import Gridap: ∇
import LinearAlgebra: tr, ⋅

const T = VectorValue{2,Float64}

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

V = TestFESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 valuetype=VectorValue{2,Float64},
 conformity=:L2)

Q = TestFESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 valuetype=Float64,
 conformity=:L2,
 constraint=:zeromean)

U = TrialFESpace(V)
P = TrialFESpace(Q)

trian = get_triangulation(model)
degree = 2*order
dΩ = LebesgueMeasure(trian,degree)

btrian = BoundaryTriangulation(model)
bdegree = 2*order
dΓ = LebesgueMeasure(btrian,bdegree)
const nb = get_normal_vector(btrian)

strian = SkeletonTriangulation(model)
sdegree = 2*order
dΛ = LebesgueMeasure(strian,sdegree)
const ns = get_normal_vector(strian)

@form a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) - ∇(q)⋅u + v⋅∇(p) )*dΩ +
  ∫( (γ/h)*v⋅u - v⋅(nb⋅∇(u)) - (nb⋅∇(v))⋅u + 2*(q*nb)⋅u )*dΓ +
  ∫(  
    (γ/h)*jump(v⊗ns)⊙jump(u⊗ns) -
      jump(v⊗ns)⊙mean(∇(u)) -
      mean(∇(v))⊙jump(u⊗ns)  +
      (γ0*h)*jump(q*ns)⋅jump(p*ns) +
      jump(q*ns)⋅mean(u) -
      mean(v)⋅jump(p*ns)
   )*dΛ

@form l((v,q)) = 
  ∫( v⋅f + q*g )*dΩ +
  ∫( (γ/h)*v⋅u - (nb⋅∇(v))⋅u + (q*nb)⋅u )*dΓ

uh, ph = solve( a((U,P),(V,Q))==l((V,Q)) )

eu = u - uh
ep = p - ph

l2(v) = v⋅v
h1(v) = v⋅v + ∇(v)⊙∇(v)

eu_l2 = sqrt(sum(∫(l2(eu))*dΩ))
eu_h1 = sqrt(sum(∫(h1(eu))*dΩ))
ep_l2 = sqrt(sum(∫(l2(ep))*dΩ))

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
