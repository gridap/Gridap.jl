module VectorValuedMultiFEOperatorsTests

# Testing Stokes problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point

using Test
using Gridap
import Gridap: ∇

const T = VectorValue{2,Float64}

# Define manufactured functions
u1(x) = VectorValue(x[1], x[2])
u2(x) = x[1] - x[2]

∇u1(x) = TensorValue(1.0,0.0,0.0,1.0)

∇(::typeof(u1)) = ∇u1

b1(x) = VectorValue(1.0,-1.0)
b2(x) = 2.0

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(2,2))

# Construct the FEspace 1
order = 2
fespace1 = FESpace(
  reffe = :QLagrangian,
  conformity = :H1,
  valuetype = T,
  model = model,
  order = order,
  diritags = "boundary")

# Construct the FEspace 2
fespace2 = FESpace(
  reffe = :PLagrangian,
  conformity = :L2,
  valuetype = Float64,
  model = model,
  order = order-1,
  constraint = :zeromean)

# Define test and trial
V1 = TestFESpace(fespace1)
V2 = TestFESpace(fespace2)
V = [V1, V2]

U1 = TrialFESpace(fespace1,u1)
U2 = TrialFESpace(fespace2)
U = [U1, U2]

# Define integration mesh and quadrature for volume
trian = Triangulation(model)
quad = CellQuadrature(trian,degree=2)

a(v,u) = 
  inner(∇(v[1]),∇(u[1])) - inner(div(v[1]),u[2]) + inner(v[2],div(u[1]))
b(v) = inner(v[1],b1) + inner(v[2],b2)
t_Ω = AffineFETerm(a,b,trian,quad)

# Define the FEOperator
op = LinearFEOperator(V,U,t_Ω)

# Solve!
uh = solve(op)

# Define exact solution and error
e1 = u1 - uh[1]

e2 = u2 - uh[2]

#writevtk(trian,"trian",cellfields=["uh2"=>uh[2],"p"=>p])

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = inner(∇(u),∇(u)) + l2(u)

# Compute errors
e1l2 = sqrt(sum( integrate(l2(e1),trian,quad) ))
e1h1 = sqrt(sum( integrate(h1(e1),trian,quad) ))

e2l2 = sqrt(sum( integrate(l2(e2),trian,quad) ))

@test e1l2 < 1.e-8
@test e1h1 < 1.e-8

@test e2l2 < 1.e-8

end # module
