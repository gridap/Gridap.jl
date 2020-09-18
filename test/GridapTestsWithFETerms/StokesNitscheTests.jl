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

# Define test and trial
Y = MultiFieldFESpace([V,Q])
# Y = [V,Q]

U = TrialFESpace(V)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
# X = [U,P]

# Define integration mesh and quadrature for volume
trian = get_triangulation(model)
quad = CellQuadrature(trian,2*order)

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,2*order)
nb = get_normal_vector(btrian)

# Weak form
function A_Ω(x,y)
  u, p = x
  v, q = y
  ∇(v)⊙∇(u) - q*(∇⋅u) - (∇⋅v)*p
end

function B_Ω(y)
  v, q = y
  v⋅f - q*g
end

function A_∂Ω(x,y)
  u, p = x
  v, q = y
  (γ/h)*v⋅u - v⋅(nb⋅∇(u)) - (nb⋅∇(v))⋅u + (p*nb)⋅v + (q*nb)⋅u
end

function B_∂Ω(y)
  v, q = y
  (γ/h)*v⋅u - (nb⋅∇(v))⋅u + (q*nb)⋅u
end

t_Ω = AffineFETerm(A_Ω,B_Ω,trian,quad)

# t_∂Ω = FESource(B_∂Ω,btrian,bquad)
t_∂Ω = AffineFETerm(A_∂Ω,B_∂Ω,btrian,bquad)

# Dummy term with 0 facets to check degenerated case
trian0 = BoundaryTriangulation(model,fill(false,Gridap.ReferenceFEs.num_facets(model)))
quad0 = CellQuadrature(trian0,2*order)
s0(x) = VectorValue(2.0*x[1],2.0)
function L0(y)
  v,q = y
  s0⋅v
end
t_0 = FESource(L0,trian0,quad0)

# Define the FEOperator
op = AffineFEOperator(X,Y,t_Ω,t_∂Ω,t_0)
# op = LinearFEOperator(Yh,Xh,t_Ω)

# Solve!
xh = solve(op)
uh, ph = xh

# Define exact solution and error
eu = u - uh

ep = p - ph

# writevtk(trian,"trian",cellfields=["uh"=>uh,"ph"=>ph, "eu"=>eu, "ep"=>ep])

# Define norms to measure the error
l2(u) = u⊙u
h1(u) = ∇(u)⊙∇(u) + l2(u)

# Compute errors
eul2 = sqrt(sum( integrate(l2(eu),quad) ))
euh1 = sqrt(sum( integrate(h1(eu),quad) ))

epl2 = sqrt(sum( integrate(l2(ep),quad) ))

@test eul2 < 1.e-8
@test euh1 < 1.e-8

@test epl2 < 1.e-8

end
