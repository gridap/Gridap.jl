module StokesNitsche

# Testing Stokes problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point

using Test
using Gridap
import Gridap: ∇

#
const T = VectorValue{2,Float64}

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

# Construct the FEspace 1
# labels = FaceLabels(model)
V = TestFESpace(
reffe=:Lagrangian,
conformity=:H1,
valuetype=VectorValue{2,Float64},
model=model,
order=order)

# Construct the FEspace 2
Q = TestFESpace(
reffe=:PLagrangian,
conformity=:L2,
valuetype=Float64,
model=model,
order=order-1,
constraint=:zeromean)

U = TrialFESpace(V)
P = TrialFESpace(Q)

# Define integration mesh and quadrature for volume
trian = get_triangulation(model)
dΩ = LebesgueMeasure(trian,2*order)

btrian = BoundaryTriangulation(model)
dΓ = LebesgueMeasure(btrian,2*order)
nb = get_normal_vector(btrian)

# Dummy term with 0 facets to check degenerated case
trian0 = BoundaryTriangulation(model,fill(false,Gridap.ReferenceFEs.num_facets(model)))
dΓ0 = LebesgueMeasure(trian0,2*order)
const s0(x) = VectorValue(2.0*x[1],2.0)

@form a((u,p),(v,q)) =
 ∫( ∇(v)⊙∇(u) - q*(∇⋅u) - (∇⋅v)*p )*dΩ +
 ∫( (γ/h)*v⋅u - v⋅(nb⋅∇(u)) - (nb⋅∇(v))⋅u + (p*nb)⋅v + (q*nb)⋅u )*dΓ

@form l((v,q)) =
 ∫( v⋅f - q*g )*dΩ +
 ∫( (γ/h)*v⋅u - (nb⋅∇(v))⋅u + (q*nb)⋅u )*dΓ +
 ∫( s0⋅v )*dΓ0

uh, ph = solve( a((U,P),(V,Q))==l((V,Q)) )

eu = u - uh
ep = p - ph

l2(v) = v⋅v
h1(v) = v⋅v + ∇(v)⊙∇(v)

eu_l2 = sqrt(sum(∫(l2(eu))*dΩ))
eu_h1 = sqrt(sum(∫(h1(eu))*dΩ))
ep_l2 = sqrt(sum(∫(l2(ep))*dΩ))

tol = 1.0e-8
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end
