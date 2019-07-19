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
res = apply(op,zh)
r2 = similar(res)
apply!(r2,op,zh)
@test res ≈ r2

cache = solve!(zh,solver,op)
@test free_dofs(zh) ≈ free_dofs(uh)

zh = zero(U)
solve!(zh,solver,op,cache)
@test free_dofs(zh) ≈ free_dofs(uh)

# With Neumann BCs

# Manufactured Neumann function
gfun(x) = 1.0

# Construct the FEspace
order = 1
diritags = [1,2,3,4,5,6,7]
fespace = ConformingFESpace(Float64,model,order,diritags)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Setup integration on Neumann boundary
neumanntags = [8,]
btrian = BoundaryTriangulation(model,neumanntags)
bquad = CellQuadrature(btrian,order=2)

# Object describing Neumann function
gfield = CellField(btrian,gfun)

# Integrand of the Neumann BC
g(v) = inner(v,gfield)

# Define weak form terms
t_Ω = AffineFETerm(a,b,trian,quad)
t_ΓN = FESource(g,btrian,bquad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define FE problem
op = LinearFEOperator(V,U,assem,t_Ω,t_ΓN)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
e = u - uh

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])

@test el2 < 1.e-8
@test eh1 < 1.e-8

# With Robin boundary conditions

# Manufactured Robin function
rfun(x) = 1.0 + ufun(x)

# Construct the FEspace
order = 1
diritags = [1,2,3,4,5,7]
fespace = ConformingFESpace(Float64,model,order,diritags)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Setup integration on Robin boundary
robintags = [6,]
rtrian = BoundaryTriangulation(model,robintags)
rquad = CellQuadrature(rtrian,order=2)

# Object describing Robin function
rfield = CellField(rtrian,rfun)

# Integrands of the Robin BC
r(v) = inner(v,rfield)
m(v,u) = inner(v,u)

# Define Robin terms
t_ΓR = AffineFETerm(m,r,rtrian,rquad)

# Dummy term that includes jumps.
# Only for testing purposes since the shape
# functions are continuous in this example.
# Thus, including this term should not change the solution
tags = [9,]
strian = SkeletonTriangulation(model,tags)
squad = CellQuadrature(strian,order=2)
j(v,u) = inner(jump(v),jump(u))
t_ΓS = LinearFETerm(j,strian,squad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define FE problem
op = LinearFEOperator(V,U,assem,t_Ω,t_ΓN,t_ΓR,t_ΓS)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
e = u - uh

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])

@test el2 < 1.e-8
@test eh1 < 1.e-8

end # module
