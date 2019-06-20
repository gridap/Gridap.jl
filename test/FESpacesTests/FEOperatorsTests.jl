module FEOperatorsTests

##
using Test
using Gridap

import Gridap: ∇
##

# Define manufactured functions
ufun(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = 0.0

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

# Define forms
a(v,u) = inner(∇(v), ∇(u))
b(v) = inner(v,bfield)

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
u = CellField(trian,ufun)
e = u - uh

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = a(u,u) + l2(u)

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])

# Further tests

@test TrialFESpace(op) === U
@test TestFESpace(op) === V

zh = zero(U)
r = apply(op,zh)
r2 = similar(r)
apply!(r2,op,zh)
@test r ≈ r2

cache = solve!(zh,solver,op)
@test free_dofs(zh) ≈ free_dofs(uh)

zh = zero(U)
solve!(zh,solver,op,cache)
@test free_dofs(zh) ≈ free_dofs(uh)

end
