module PoissonTests

using Test
using Gridap
import Gridap: ∇
using LinearAlgebra

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)
order = 2

const h = (domain[2]-domain[1]) / partition[1]
const γ = 10

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[7,8])
add_tag_from_tags!(labels,"nitsche",6)

trian = get_triangulation(model)
degree = order
dΩ = LebesgueMeasure(trian,degree)

ntrian = BoundaryTriangulation(model,labels,"neumann")
ndegree = order
dΓn = LebesgueMeasure(ntrian,ndegree)
const nn = get_normal_vector(ntrian)

dtrian = BoundaryTriangulation(model,labels,"nitsche")
ddegree = order
dΓd = LebesgueMeasure(dtrian,ddegree)
const dn = get_normal_vector(dtrian)

# Using automatic differentiation
u_scal(x) = x[1]^2 + x[2]
f_scal(x) = - Δ(u_scal)(x)

#u_scal(x) = x[1]^2 + x[2]
#∇u_scal(x) = VectorValue( 2*x[1], one(x[2]) )
#Δu_scal(x) = 2
#f_scal(x) = - Δu_scal(x)
#∇(::typeof(u_scal)) = ∇u_scal

scalar_data = Dict{Symbol,Any}()
scalar_data[:valuetype] = Float64
scalar_data[:u] = u_scal
scalar_data[:f] = f_scal

# Using automatic differentiation
u_vec(x) = VectorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2 )
f_vec(x) = - Δ(u_vec)(x)

#u_vec(x) = VectorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2 )
#∇u_vec(x) = TensorValue( 2*x[1], one(x[2]), 4*one(x[1]), - 2*x[2] )
#Δu_vec(x) = VectorValue( 2, -2 )
#f_vec(x) = - Δu_vec(x)
#∇(::typeof(u_vec)) = ∇u_vec

vector_data = Dict{Symbol,Any}()
vector_data[:valuetype] = VectorValue{2,Float64}
vector_data[:u] = u_vec
vector_data[:f] = f_vec

for data in [ scalar_data, vector_data]

  T = data[:valuetype]
  u = data[:u]
  f = data[:f]

  V = TestFESpace(
   model=model,
   order=order,
   reffe=:Lagrangian,
   labels=labels,
   valuetype=T,
   dirichlet_tags="dirichlet")

  U = TrialFESpace(V,u)

  uh = interpolate(u, U)

  @form a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ + ∫( (γ/h)*v⊙u  - v⊙(dn⋅∇(u)) - (dn⋅∇(v))⊙u )*dΓd

  @form l(v) = ∫( v⊙f )*dΩ + ∫( v⊙(nn⋅∇(uh)) )*dΓn + ∫( (γ/h)*v⊙uh - (dn⋅∇(v))⊙u )*dΓd

  uh = solve( a(U,V) == l(V) )

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

end

end # module
