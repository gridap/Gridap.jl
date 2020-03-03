module StokesNitsche

# Testing Stokes problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point

using Test
using Gridap
import Gridap: ∇
import LinearAlgebra: tr

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
  inner(∇(v), ∇(u)) - inner(q,divergence(u)) - inner(divergence(v), p)
end

function B_Ω(y)
  v, q = y
  inner(v,f) - inner(q, g)
end

function A_∂Ω(x,y)
  u, p = x
  v, q = y
  # (γ/h) * inner(v,u) - inner(outer(nb,v), ∇(u)) - inner(∇(v), outer(nb,u)) + inner(v, p*nb) + inner(q*nb,u)
  (γ/h)*v*u - v*(nb*∇(u)) - (nb*∇(v))*u + (p*nb)*v + (q*nb)*u
end

function B_∂Ω(y)
  v, q = y
  # + (γ/h) * inner(v,ud)  - inner(∇(v), outer(nb,ud_cf)) + inner(q*nb,ud)
  (γ/h)*v*u - (nb*∇(v))*u + (q*nb)*u
end

t_Ω = AffineFETerm(A_Ω,B_Ω,trian,quad)

# t_∂Ω = FESource(B_∂Ω,btrian,bquad)
t_∂Ω = AffineFETerm(A_∂Ω,B_∂Ω,btrian,bquad)

# Define the FEOperator
op = AffineFEOperator(X,Y,t_Ω,t_∂Ω)
# op = LinearFEOperator(Yh,Xh,t_Ω)

# Solve!
xh = solve(op)
uh, ph = xh

# Define exact solution and error
eu = u - uh

ep = p - ph

# writevtk(trian,"trian",cellfields=["uh"=>uh,"ph"=>ph, "eu"=>eu, "ep"=>ep])

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = inner(∇(u),∇(u)) + l2(u)

# Compute errors
eul2 = sqrt(sum( integrate(l2(eu),trian,quad) ))
euh1 = sqrt(sum( integrate(h1(eu),trian,quad) ))

epl2 = sqrt(sum( integrate(l2(ep),trian,quad) ))

@test eul2 < 1.e-8
@test euh1 < 1.e-8

@test epl2 < 1.e-8

end
