module VectorValuedMultiFEOperatorsTests

# Testing Stokes problem with Dirichlet BCs for the velocity
# and pressure prescribed in a single point

using Test
using Gridap
import Gridap: ∇

const T = VectorValue{2,Float64}

# Define manufactured functions
u1fun(x) = VectorValue(x[1], x[2])
u2fun(x) = x[1] - x[2]

u1fun_grad(x) = TensorValue(1.0,0.0,0.0,1.0)

∇(::typeof(u1fun)) = u1fun_grad

b1fun(x) = VectorValue(1.0,-1.0)
b2fun(x) = 2.0

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace 1
order = 2
diritag = "boundary"
fespace1 = CLagrangianFESpace(T,model,order,diritag)

# Construct the FEspace 2
diritag = 1
fespace2 = CLagrangianFESpace(Float64,model,order-1,diritag)

# Define test and trial
V1 = TestFESpace(fespace1)
V2 = TestFESpace(fespace2)
V = [V1, V2]

U1 = TrialFESpace(fespace1,u1fun)
U2 = TrialFESpace(fespace2,u2fun)
U = [U1, U2]

# Define integration mesh and quadrature for volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

a(v,u) = 
  inner(∇(v[1]),∇(u[1])) - inner(div(v[1]),u[2]) + inner(v[2],div(u[1]))
b(v) = inner(v[1],b1fun) + inner(v[2],b2fun)
t_Ω = AffineFETerm(a,b,trian,quad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(V,U,assem,t_Ω)

# Define the FESolver
ls = LUSolver()
solver = LinearFESolver(ls)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
u1 = CellField(trian,u1fun)
e1 = u1 - uh[1]

u2 = CellField(trian,u2fun)
e2 = u2 - uh[2]

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
