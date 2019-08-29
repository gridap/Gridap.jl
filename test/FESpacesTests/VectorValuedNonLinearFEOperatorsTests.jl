module VectorValuedNonLinearFEOperatorsTests

using Test
using Gridap
import Gridap: ∇

# Material parameters
const λ = 30.0
const μ = 40.0

# Constitutive law (St. Venant–Kirchhoff Material)
S(E) = λ*tr(E)*one(E) + 2*μ*E

# Green strain
E(F) = 0.5*( F'*F - one(F) )
dE(F,dF) = 0.5*( dF'*F + F'*dF)

# Operations to be performed at integration points
@law σ(x,∇u) = ∇u * S(E(∇u))
@law dσ(x,∇du,∇u) = ∇du * S(E(∇u)) + ∇u * S(dE(∇u,∇du))

# Define manufactured functions
ufun(x) = VectorValue(x[1] + x[2],x[1])
ufun_grad(x) = TensorValue(1.0,1.0,1.0,0.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = VectorValue(0.0,0.0)

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 2
diritag = "boundary"
T = VectorValue{2,Float64}
fespace = ConformingFESpace(T,model,order,diritag)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Setup integration on the volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

# Transform shape function according to the model

# Terms in the volume
bfield = CellField(trian,bfun)
res(u,v) = inner( ∇(v), σ(∇(u)) ) - inner(v,bfield)
jac(u,v,du) = inner(∇(v), dσ(∇(du),∇(u)) )
t_Ω = NonLinearFETerm(res,jac,trian,quad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = NonLinearFEOperator(V,U,assem,t_Ω)

# Define the FESolver
ls = LUSolver()
tol = 1.e-10
maxiters = 20
nls = NewtonRaphsonSolver(ls,tol,maxiters)
solver = NonLinearFESolver(nls)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
u = CellField(trian,ufun)
e = u - uh

# Define norms to measure the error
l2(u) = inner(u,u)
sh1(u) = inner(∇(u),∇(u))
h1(u) = sh1(u) + l2(u)

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8


# Further tests
op = NonLinearFEOperator(V,U,t_Ω)
op = NonLinearFEOperator(V,U,t_Ω,t_Ω)

end # module
