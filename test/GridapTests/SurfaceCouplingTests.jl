module SurfaceCouplingTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.FESpaces
import Gridap: ∇
using LinearAlgebra: tr, ⋅

# Analytical functions

u(x) = VectorValue( x[1]^2 + 2*x[2]^2, -x[1]^2 )
∇u(x) = TensorValue( 2*x[1], 4*x[2], -2*x[1], zero(x[1]) )
Δu(x) = VectorValue( 6, -2 )

p(x) = x[1] + 3*x[2]
∇p(x) = VectorValue(1,3)

s(x) = -Δu(x)
f(x) = -Δu(x) + ∇p(x)
g(x) = tr(∇u(x))

∇(::typeof(u)) = ∇u
∇(::typeof(p)) = ∇p

# Geometry + Integration

n = 20
mesh = (n,n)
domain = 2 .* (0,1,0,1) .- 1
order = 1
model = CartesianDiscreteModel(domain, mesh)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[6,7,8])

trian = Triangulation(model)

const R = 0.4

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

cell_to_coods = get_cell_coordinates(trian)
cell_to_is_solid = collect1d(apply(is_in,cell_to_coods))
cell_to_is_fluid = Vector{Bool}(.! cell_to_is_solid)

trian_solid = Triangulation(model, cell_to_is_solid)
trian_fluid = Triangulation(model, cell_to_is_fluid)

order = 2

degree = 2*order
dΩ = LebesgueMeasure(trian,degree)
dΩ_solid = LebesgueMeasure(trian_solid,degree)
dΩ_fluid = LebesgueMeasure(trian_fluid,degree)

btrian = BoundaryTriangulation(model,labels,"neumann")
bdegree = 2*order
dΛ = LebesgueMeasure(btrian,bdegree)
n = get_normal_vector(btrian)

# This returns a SkeletonTriangulation whose normal vector
# goes outwards to the fluid domain.
trian_Γ = InterfaceTriangulation(model,cell_to_is_fluid)
n_Γ = get_normal_vector(trian_Γ)
dΓ = LebesgueMeasure(trian_Γ,bdegree)

# FESpaces

V = TestFESpace(
  model=model,
  valuetype=VectorValue{2,Float64},
  reffe=:QLagrangian,
  order=order,
  conformity =:H1,
  dirichlet_tags="dirichlet")

Q = TestFESpace(
  triangulation=trian_fluid,
  valuetype=Float64,
  order=order-1,
  reffe=:PLagrangian,
  conformity=:L2)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

# FE problem

@form a((u,p),(v,q)) =
  ∫( inner(∇(v),∇(u)) )*dΩ_solid +
  ∫( inner(∇(v),∇(u)) - (∇⋅v)*p + q*(∇⋅u) )*dΩ_fluid

@form l((v,q)) =
  ∫( v⋅s )*dΩ_solid +
  ∫( v⋅f + q*g )*dΩ_fluid +
  ∫( v⋅(n⋅∇u) - (n⋅v)*p )*dΛ +
  ∫( - mean(n_Γ⋅v)*p )*dΓ

uh, ph = solve( a((U,P),(V,Q))==l((V,Q)) )

eu = u - uh
ep = p - ph

# Visualization

d = mktempdir()
writevtk(trian,joinpath(d,"trian"),cellfields=["ep"=>ep,"eu"=>eu,"p"=>p])
writevtk(trian_fluid,joinpath(d,"trian_fluid"),cellfields=["ep"=>ep,"p"=>p,"ph"=>ph])

# Errors

l2(v) = v⋅v
h1(v) = v⋅v + ∇(v)⊙∇(v)

eu_l2 = sqrt(sum(∫(l2(eu))*dΩ))
eu_h1 = sqrt(sum(∫(h1(eu))*dΩ))
ep_l2 = sqrt(sum(∫(l2(ep))*dΩ_fluid))

tol = 1.0e-8
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
