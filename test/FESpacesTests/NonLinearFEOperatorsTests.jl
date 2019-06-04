module NonLinearFEOperatorsTests

using Gridap.FESpaces
using Gridap.Assemblers
using Gridap.FEOperators
using Gridap.LinearSolvers
using Gridap.NonLinearSolvers

using Test
using Gridap
using Gridap.CellMaps
using Gridap.Geometry
using Gridap.Geometry.Cartesian
using Gridap.FieldValues
using Gridap.CellQuadratures
using Gridap.CellIntegration
using Gridap.Vtkio

import Gridap: gradient

# Define manufactured functions
ufun(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0)
gradient(::typeof(ufun)) = ufun_grad
bfun(x) = -(3.0*x[1]+x[2]+1.0)
νfun(x,u) = (u+1.0)*x[1]

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

# Define cell field describing the source term
bfield = CellField(trian,bfun)

# Define a solution dependent material parameter
ν(u) = CellField(trian,νfun,u)

# Define residual and jacobian
a(u,v,du) = inner( ∇(v), ν(u)*∇(du))
res(u,v) = a(u,v,u) - inner(v,bfield)
jac(u,v,du) = a(u,v,du)  # + inner(v,ν(du)*grad(u))

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
u = CellField(trian,ufun)
e = u - uh

# Define norms to measure the error
l2(u) = inner(u,u)
sh1(u) = inner(∇(u),∇(u))
h1(u) = sh1(u) + l2(u)

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
ul2 = sqrt(sum( integrate(l2(u),trian,quad) ))
uh1 = sqrt(sum( integrate(h1(u),trian,quad) ))

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])

# Further tests
@test TrialFESpace(op) === U
@test TestFESpace(op) === V

zh = zero(U)
cache = solve!(zh,solver,op)
@test free_dofs(zh) ≈ free_dofs(uh)

zh = zero(U)
solve!(zh,solver,op,cache)
@test free_dofs(zh) ≈ free_dofs(uh)

end # module NonLinearFEOperatorsTests
