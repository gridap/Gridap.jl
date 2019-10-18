module NewFESpacesTests
##
using Gridap.DOFBases: LagrangianDOFBasis

using Gridap, Test

using Gridap.Helpers

using Base.Cartesian

import Gridap: ∇
# 1D reffe

F = Gridap.NewConformingFESpaces

# Construct the discrete model
n = 2
model = CartesianDiscreteModel(partition=(n,))

diri = [1]
p = Polytope(1,)
order = 2

trian = Triangulation(model)
reffe = LagrangianRefFE(Float64,p,order)

graph = GridGraph(model)
labels = FaceLabels(model)
mydiri_tags = Gridap.BoundaryGrids._setup_tags(labels,diri)

fespace = ConformingFESpace(reffe,trian,graph,labels,mydiri_tags)
# Define manufactured functions
ufun(x) = x[1]*(x[1]-2.0)
ufun_grad(x) = VectorValue(2*x[1]-2.0,)
bfun(x) = -2.0
∇(::typeof(ufun)) = ufun_grad

V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

##
# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,degree=order*2)

# Define forms
a(v,u) = inner(∇(v), ∇(u))
b(v) = inner(v,bfun)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(a,b,V,U,assem,trian,quad)

# Define the FESolver
ls = BackslashSolver()
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

# writevtk(trian,"onedres",cellfields=["uh"=>uh,"e"=>e])
##
end # module
