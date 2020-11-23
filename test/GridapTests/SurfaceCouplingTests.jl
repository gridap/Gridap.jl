module SurfaceCouplingTests

using Test
using Gridap
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

Ω = Triangulation(model)

const R = 0.4

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

cell_to_coods = get_cell_coordinates(Ω)
cell_to_is_solid = lazy_map(is_in,cell_to_coods)
cell_to_is_fluid = lazy_map(!,cell_to_is_solid)

model_solid = DiscreteModel(model,cell_to_is_solid)
model_fluid = DiscreteModel(model,cell_to_is_fluid)

Ωs = Triangulation(model_solid)
Ωf = Triangulation(model_fluid)
Λ = BoundaryTriangulation(model,labels,tags="neumann")
Γ = InterfaceTriangulation(model_fluid,model_solid)

n_Λ = get_normal_vector(Λ)
n_Γ = get_normal_vector(Γ)

order = 2
degree = 2*order

dΩ = LebesgueMeasure(Ω,degree)
dΩs = LebesgueMeasure(Ωs,degree)
dΩf = LebesgueMeasure(Ωf,degree)
dΛ = LebesgueMeasure(Λ,degree)
dΓ = LebesgueMeasure(Γ,degree)

# FE Spaces

reffe_u = ReferenceFE(:Lagrangian,VectorValue{2,Float64},order)
reffe_p = ReferenceFE(:Lagrangian,Float64,order-1,space=:P)

V = TestFESpace(model,reffe_u,conformity=:H1,labels=labels,dirichlet_tags="dirichlet")
Q = TestFESpace(model_fluid,reffe_p,conformity=:L2)
U = TrialFESpace(V,u)
P = Q

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

#uh, ph = FEFunction(X,rand(num_free_dofs(X)))
#vh, qh = FEFunction(Y,rand(num_free_dofs(Y)))
#writevtk(Ω,"trian",cellfields=["uh"=>uh,"ph"=>ph,"vh"=>vh,"qh"=>qh])

# Weak form

a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) )*dΩs +
  ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )*dΩf

l((v,q)) =
  ∫( v⋅s )*dΩs +
  ∫( v⋅f + q*g )*dΩf +
  ∫( v⋅(n_Λ⋅∇u) - (n_Λ⋅v)*p )*dΛ +
  ∫( - (n_Γ.⁺⋅v.⁺)*p )*dΓ

# FE problem

op = AffineFEOperator(a,l,X,Y)
uh, ph = solve(op)

# Visualization

eu = u - uh
ep = p - ph

#writevtk(Ω,"trian",cellfields=["uh"=>uh,"ph"=>ph,"eu"=>eu,"ep"=>ep])
#writevtk(Ωf,"trian_fluid",cellfields=["uh"=>uh,"ph"=>ph,"eu"=>eu,"ep"=>ep])

# Errors

eu_l2 = sqrt(sum(∫( eu⋅eu )*dΩ))
eu_h1 = sqrt(sum(∫( eu⋅eu + ∇(eu)⊙∇(eu) )*dΩ))
ep_l2 = sqrt(sum(∫( ep*ep )*dΩf))

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
