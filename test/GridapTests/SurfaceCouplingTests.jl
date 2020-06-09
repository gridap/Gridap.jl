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
quad = CellQuadrature(trian,degree)
quad_solid = CellQuadrature(trian_solid,degree)
quad_fluid = CellQuadrature(trian_fluid,degree)

btrian = BoundaryTriangulation(model,labels,"neumann")
bdegree = 2*order
bquad = CellQuadrature(btrian,bdegree)
n = get_normal_vector(btrian)

# This returns a SkeletonTriangulation whose normal vector
# goes outwards to the fluid domain.
trian_Γ = InterfaceTriangulation(model,cell_to_is_fluid)
n_Γ = get_normal_vector(trian_Γ)
quad_Γ = CellQuadrature(trian_Γ,bdegree)

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

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

# FE problem

function a_solid(x,y)
  u,p = x
  v,q = y
  inner(∇(v),∇(u))
end

function l_solid(y)
  v,q = y
  v⋅s
end

function a_fluid(x,y)
  u,p = x
  v,q = y
  inner(∇(v),∇(u)) - (∇⋅v)*p + q*(∇⋅u)
end

function l_fluid(y)
  v,q = y
  v⋅f + q*g
end

function l_Γn_fluid(y)
  v,q = y
  v⋅(n⋅∇u) - (n⋅v)*p
end

# Pressure drop at the interface
function l_Γ(y)
  v,q = y
  - mean(n_Γ⋅v)*p
end

t_Ω_solid = AffineFETerm(a_solid,l_solid,trian_solid,quad_solid)
t_Ω_fluid = AffineFETerm(a_fluid,l_fluid,trian_fluid,quad_fluid)
t_Γn_fluid = FESource(l_Γn_fluid,btrian,bquad)
t_Γ = FESource(l_Γ,trian_Γ,quad_Γ)

op = AffineFEOperator(X,Y,t_Ω_solid,t_Ω_fluid,t_Γn_fluid,t_Γ)
uh, ph = solve(op)

# Visualization

ph_fluid = restrict(ph, trian_fluid)

eu = u - uh
ep = p - ph
ep_fluid = p - ph_fluid

#writevtk(trian_fluid,"trian_fluid",cellfields=["ph"=>ph_fluid, "ep"=>ep_fluid])
#
#writevtk(trian,"trian", cellfields=["uh" => uh, "ph"=> ph, "eu"=>eu, "ep"=>ep])

# Errors

l2(v) = v⋅v
h1(v) = v⋅v + inner(∇(v),∇(v))

eu_l2 = sqrt(sum(integrate(l2(eu),trian,quad)))
eu_h1 = sqrt(sum(integrate(h1(eu),trian,quad)))
ep_l2 = sqrt(sum(integrate(l2(ep_fluid),trian_fluid,quad_fluid)))

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
