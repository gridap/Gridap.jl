module PoissonDGTests

using Test
using Gridap
import Gridap: ∇
using LinearAlgebra

#domain = (0,1,0,1)
#partition = (4,4)
#model = CartesianDiscreteModel(domain,partition)
#const h = (domain[2]-domain[1]) / partition[1]

using Gridap.Geometry: DiscreteModelMock
model = DiscreteModelMock()
const h = 1

order = 2
const γ = 10

trian = get_triangulation(model)
degree = 2*order
dΩ = LebesgueMeasure(trian,degree)

btrian = BoundaryTriangulation(model)
bdegree = 2*order
dΓb = LebesgueMeasure(btrian,bdegree)
const bn = get_normal_vector(btrian)

strian = SkeletonTriangulation(model)
sdegree = 2*order
dΓs = LebesgueMeasure(strian,sdegree)
const sn = get_normal_vector(strian)

u_scal(x) = x[1]^2 + x[2]
∇u_scal(x) = VectorValue( 2*x[1], one(x[2]) )
Δu_scal(x) = 2
f_scal(x) = - Δu_scal(x)
∇(::typeof(u_scal)) = ∇u_scal

T = Float64
u = u_scal
f = f_scal

V = TestFESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 conformity=:L2,
 valuetype=T)

U = TrialFESpace(V,u)

@form a(u,v) =
  ∫( ∇(v)⋅∇(u) )*dΩ +
  ∫( (γ/h)*v*u  - v*(bn⋅∇(u)) - (bn⋅∇(v))*u )*dΓb +
  ∫( (γ/h)*jump(v*sn)⋅jump(u*sn) - jump(v*sn)⋅mean(∇(u)) - mean(∇(v))⋅jump(u*sn) )*dΓs

@form l(v) =
  ∫( v*f )*dΩ +
  ∫( (γ/h)*v*u - (bn⋅∇(v))*u )*dΓb

uh = solve( a(U,V) == l(V) )

e = u - uh

l2(u) = u⊙u
sh1(u) = ∇(u)⊙∇(u)
h1(u) = sh1(u) + l2(u)

el2 = sqrt(sum( ∫( l2(e) )*dΩ ))
eh1 = sqrt(sum( ∫( h1(e) )*dΩ ))
ul2 = sqrt(sum( ∫( l2(uh) )*dΩ ))
uh1 = sqrt(sum( ∫( h1(uh) )*dΩ ))

d = mktempdir()

writevtk(trian,joinpath(d,"trian"),order=order,
  cellfields=["uh"=>uh,"e"=>e,"u"=>u],
  celldata=["el2"=>∫( l2(e) )*dΩ])

writevtk(strian,joinpath(d,"strian"),order=order,cellfields=[
  "jump_e"=>jump(e),
  "jump_∇e"=>jump(∇(e)),
  "jump_n∇e"=>jump(sn⋅∇(e)),
  "jump_∇uh"=>jump(∇(uh)),
  "jump_n∇uh"=>jump(sn⋅∇(uh)),
  "jump_uh"=>jump(uh)])

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
