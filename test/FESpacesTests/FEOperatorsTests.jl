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

include("FEOperatorsTestsMixin.jl")

# TODO Following test not yet working for triangles
# since it involves jumps and they do not work for non-oriented facets

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

# Further tests
op = LinearFEOperator(V,U,t_Ω,t_ΓN,t_ΓR,t_ΓS)
op = LinearFEOperator(V,U,t_Ω)

# For triangles

model = simplexify(model)

include("FEOperatorsTestsMixin.jl")

end # module
