module StokesNitsche

# Testing Stokes problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point

using Test
using Gridap
import Gridap: ∇

u(x) = VectorValue(x[1]*x[1], x[2])
∇u(x) = TensorValue(2*x[1],0.0,0.0,1.0)
Δu(x) = VectorValue(2.0,0.0)

p(x) = x[1] - x[2]
∇p(x) = VectorValue(1.0,-1.0)

∇(::typeof(u)) = ∇u

n = 2
order = 2

# function run_test(n,order)

f(x) = - Δu(x) + ∇p(x)
g(x) = tr(∇u(x)) # i.e. div(u(x))
ud(x) = u(x)

# Discrete model
L = 1.0
limits = (0.0, L, 0.0, L)
ncellx = n
model = CartesianDiscreteModel(limits, (ncellx,ncellx))

h = L / ncellx

γ = order*(order+1)
γ0 = 1.0/10.0

s0(x) = VectorValue(2.0*x[1],2.0)

reffe_u = ReferenceFE(:Lagrangian,VectorValue{2,Float64},order)
reffe_p = ReferenceFE(:Lagrangian,Float64,order-1,space=:P)

U = FESpace(model,reffe_u,conformity=:H1)
P = FESpace(model,reffe_p,conformity=:L2,constraint=:zeromean)

X = MultiFieldFESpace([U,P])

# Define integration mesh and quadrature for volume
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
n_Γ = get_normal_vector(Γ)
# Dummy term with 0 facets to check degenerated case
Γ0 = BoundaryTriangulation(model,fill(false,num_facets(model)))

# Lebesgue measures
degree = 2*order
dΩ = LebesgueMeasure(Ω,degree)
dΓ = LebesgueMeasure(Γ,degree)
dΓ0 = LebesgueMeasure(Γ0,degree)

# Weak form

a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) - q*(∇⋅u) - (∇⋅v)*p )*dΩ +
  ∫( (γ/h)*v⋅u - v⋅(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))⋅u + (p*n_Γ)⋅v + (q*n_Γ)⋅u )*dΓ

l((v,q)) =
  ∫( v⋅f - q*g )*dΩ +
  ∫( (γ/h)*v⋅u - (n_Γ⋅∇(v))⋅u + (q*n_Γ)⋅u )*dΓ + ∫( s0⋅v )*dΓ0

# Define the FEOperator
op = AffineFEOperator(a,l,X,X)

# Solve!
xh = solve(op)
uh, ph = xh

# Define exact solution and error
eu = u - uh

ep = p - ph

# writevtk(trian,"trian",cellfields=["uh"=>uh,"ph"=>ph, "eu"=>eu, "ep"=>ep])

# Define norms to measure the error
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
