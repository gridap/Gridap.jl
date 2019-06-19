module MultiFEOperatorsTests

using Test
using Gridap
import Gridap: ∇

# Define manufactured functions
u1fun(x) = x[1] + x[2]
u2fun(x) = x[1] - x[2]

u1fun_grad(x) = VectorValue(1.0,1.0)
u2fun_grad(x) = VectorValue(1.0,-1.0)

∇(::typeof(u1fun)) = u1fun_grad
∇(::typeof(u2fun)) = u2fun_grad

b1fun(x) = u2fun(x)
b2fun(x) = 0.0

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

# Define test and trial
V1 = TestFESpace(fespace)
V2 = V1
V = [V1, V2]

U1 = TrialFESpace(fespace,u1fun)
U2 = TrialFESpace(fespace,u2fun)
U = [U1, U2]

# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

# Define cell field describing the source term
b1field = CellField(trian,b1fun)
b2field = CellField(trian,b2fun)

# Define forms
a(v,u) = inner(∇(v[1]),∇(u[1])) + inner(v[1],u[2]) + inner(∇(v[2]),∇(u[2]))

b(v) = inner(v[1],b1field) + inner(v[2],b2field)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(a,b,V,U,assem,trian,quad)

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
e2h1 = sqrt(sum( integrate(h1(e2),trian,quad) ))

@test e1l2 < 1.e-8
@test e1h1 < 1.e-8

@test e2l2 < 1.e-8
@test e2h1 < 1.e-8

end
