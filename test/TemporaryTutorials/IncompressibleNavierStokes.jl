module IncompressibleNavierStokes

using Test
using Gridap
# The gradient method, also represented with ∇, must be extended
# for the functions being defined by the user
import Gridap: ∇

# Define manufactured functions
ufun(x) = VectorValue(x[1]^2,x[1]-x[2])
pfun(x) = 2*x[1]

# Gradient of manufactured solutions
ufun_grad(x) = TensorValue(2*x[1],1.0,0.0,-1.0)
pfun_grad(x) = VectorValue(2.0)
# Tensors are column-major, as in Julia

# Here is where we define gradient for our previously defined functions
∇(::typeof(ufun)) = ufun_grad
∇(::typeof(pfun)) = pfun_grad

# Body force
bfun_u(x) = VectorValue(0.0,0.0)
bfun_p(x) = 0.0

# Construct the discrete model (geometry)
# The domain is [0,1]^2 and the partition is 4x4
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 2
# All vertices and edges in the square to have Dirichlet data
diritags = [1,2,3,4,5,6,7,8]
# Only enforce the x1-component in the upper edge
# We must provide for every geometrical entity in diritags which
# are the components to be fixed in a tuple of Booleans
dirimasks = Tuple{Bool,Bool}[]
for i in 1:7
  push!(dirimasks,(true,true))
end
push!(dirimasks,(true,false))
# Vector field for velocity
T = VectorValue{2,Float64}
u_fesp = CLagrangianFESpace(T,model,order,diritags,dirimasks)

# Now we create the pressure space, one order less and discontinuous
# @sbadia : I would eliminate the need of diritags and dirimasks, e.g., dirimasks
# for scalar fields do not has sense
diritags = Int[]
dirimasks = Tuple{Bool,Bool}[]
T = Float64
p_fesp = DLagrangianFESpace(T,model,order-1,diritags, dirimasks)

# Define test and trial
# The TestFESpace method takes a FESpace and enforces homogeneous bc's on
# Dirichlet boundary
V = TestFESpace(u_fesp)
Q = TestFESpace(p_fesp)
Y = [V, Q]

U = TrialFESpace(u_fesp,ufun)
# I think we do not need to provide a function if not needed
P = TrialFESpace(p_fesp,pfun)
X = [U, P]

# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

# Define cell field describing the source term
b_u = CellField(trian,bfun_u)
b_p = CellField(trian,bfun_p)
b_tot = [b_u, b_p]
# Define a solution dependent material parameter
# ν(u1) = CellField(trian,νfun,u1)
# dν(du1) = CellBasis(trian,dνfun,du1)

# Define residual and jacobian
ν = 1.0
# + inner(v[1],u[1]*∇u[1])
# inner(∇(v[1]),ν*∇(b_u))
a(u,v) = inner(∇(v[1]),∇(u[1])) - inner(tr(∇(v[1])),u[2]) + inner(v[2],tr(∇(u[1])))
da(u,v,du) = a(du,v)
b(v) = inner(v[1],b_u) + inner(v[2],b_p)

res(u,v) = a(u,v) - b(v)
jac(u,v,du) = da(u,v,du)

# Define Assembler
assem = SparseMatrixAssembler(Y,X)

# Define the FEOperator
op = NonLinearFEOperator(res,jac,Y,X,assem,trian,quad)

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

end



# D = 2
# dirimasks = Vector{NTuple{D,Bool}}(undef,length(diritags))
# flag = ntuple(i -> true, D)
# fill!(dirimasks,flag)
