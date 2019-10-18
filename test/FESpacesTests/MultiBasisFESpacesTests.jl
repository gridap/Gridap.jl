module MultiBasisFESpaces
##
using Gridap.DOFBases: LagrangianDOFBasis

using Gridap, Test

using Gridap.Helpers

using Base.Cartesian

import Gridap: ∇
# 1D reffe

include("ScaledNodesLagrangianRefFEs.jl")
include("MultiBasisFESpaces.jl")
F = Gridap.ConformingFESpaces

D = 1
n = 5
ℓ = 1.0
β = 4.0
h = 1.0 / n
h0 = 1.0 / (n-1+β)
r = h/h0
e = h/(ℓ-h0*(n-1))
ℓ0 = n*h
geomap = ( x -> ( x[1] <= (n-1)*h ) ? x[1]/r : (n-1)*h0+(x[1]-(n-1)*h)/e )
geomap = ( x -> x[1] <= (ℓ-h) ? x[1] : (x[1]-ℓ+h)*β+ℓ-h)
newl = geomap(1.0)
typeof(newl)
# geomap = ( x -> x )
model = CartesianDiscreteModel(partition=(n,),map=geomap)
model.ugrid
# Construct the discrete model

# model = CartesianDiscreteModel(partition=(2,))
diri = [1]
p = Polytope(1,)
order = 2

trian = Triangulation(model)
reffe = LagrangianRefFE(Float64,p,order)
sreffe = ScaledNodesLagrangianRefFE(Float64,p,order,1.0/β)
values = LagrangianRefFE{1,Float64}[reffe,sreffe]
pointer = [ i == ncells(trian) ? 2 : 1 for i in 1:ncells(trian)]
reffes = CompressedCellValue(values,pointer)
typeof(reffes)

# grid = Grid(model,D)
# trian = Triangulation(grid)
# orders = fill(order,D)
graph = GridGraph(model)
labels = FaceLabels(model)
mydiri_tags = Gridap.BoundaryGrids._setup_tags(labels,diri)

# fespace = ConformingFESpace(reffe,trian,graph,labels,[4])
fespace = MultiShapeConformingFESpace(reffes,trian,graph,labels,mydiri_tags)
# typeof(reffes)
# isa(reffes,IndexCellValue{LagrangianRefFE{1,Float64},1})

# Define manufactured functions
newl = 1.6
# @santiagobadia : I do not understand what is going on
# Isn't it a closure?
# ufun(x) = x[1]*(x[1]-2.0*newl)
# ufun_grad(x) = VectorValue(2*x[1]-newl,)
ufun(x) = x[1]*(x[1]-2.0*1.6)
ufun_grad(x) = VectorValue(2*x[1]-2.0*1.6,)
bfun(x) = -2.0
# ufun(x) = x[1]
# ufun_grad(x) = VectorValue(1.0,)
# bfun(x) = 0.0
∇(::typeof(ufun)) = ufun_grad

V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)
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
uh.free_dofs
uh.diri_dofs


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
