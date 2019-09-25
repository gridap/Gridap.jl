module MultiNonLinearFEOperatorsTests

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

b1fun(x) = u2fun(x) -(3.0*x[1]+x[2]+1.0)
b2fun(x) = 0.0

@law ν(x,u1) = (u1+1.0)*x[1]
@law dν(x,du1) = du1*x[1]

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 2
diritag = "boundary"
fespace = H1ConformingFESpace(Float64,model,order,diritag)

# Define test and trial
V1 = TestFESpace(fespace)
V2 = V1
V = [V1, V2]

U1 = TrialFESpace(fespace,u1fun)
U2 = TrialFESpace(fespace,u2fun)
U = [U1, U2]

# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,degree=3*order-1)

# Define cell field describing the source term
b1field = CellField(trian,b1fun)
b2field = CellField(trian,b2fun)

# Define residual and jacobian
a(u,v,du) = inner(∇(v[1]),ν(u[1])*∇(du[1])) + inner(v[1],du[2]) + inner(∇(v[2]),∇(du[2]))
b(v) = inner(v[1],b1field) + inner(v[2],b2field)

res(u,v) = a(u,v,u) - b(v)
jac(u,v,du) = a(u,v,du) + inner(∇(v[1]),dν(du[1])*∇(u[1]))

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = NonLinearFEOperator(res,jac,V,U,assem,trian,quad)

# Define the FESolver
ls = LUSolver()
tol = 1.e-10
maxiters = 20
nls = NewtonRaphsonSolver(ls,tol,maxiters)
solver = NonLinearFESolver(nls)

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

# Further tests

t_Ω = NonLinearFETerm(res,jac,trian,quad)
op = NonLinearFEOperator(V,U,t_Ω)

uh = solve(solver,op)


end
