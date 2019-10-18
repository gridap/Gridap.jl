module NewFEFESpacesTests
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

fesp = NewConformingFESpace(reffe,trian,graph,labels,mydiri_tags)
##
# fespace = ConformingFESpace(reffe,trian,graph,labels,[4])
# fespace = MultiShapeConformingFESpace(reffes,trian,graph,labels,mydiri_tags)
# typeof(reffes)
# isa(reffes,IndexCellValue{LagrangianRefFE{1,Float64},1})

# Define manufactured functions
ufun(x) = x[1]*(x[1]-2.0)
ufun_grad(x) = VectorValue(2*x[1]-2.0,)
bfun(x) = -2.0
# ufun(x) = x[1]
# ufun_grad(x) = VectorValue(1.0,)
# bfun(x) = 0.0
∇(::typeof(ufun)) = ufun_grad

fun = ufun

T = Float64

reffe = fesp.reffe
dofb = broadcast(dofbasis,reffe)
trian = fesp.triangulation
phi = CellGeomap(trian)
uphys = fun ∘ phi
celldofs = fesp.cell_eqclass
nfdofs = fesp.dim_to_nface_eqclass
E = dof_type(T)
free_dofs = zeros(E, num_free_dofs(fesp))
diri_dofs = zeros(E, num_diri_dofs(fesp))
aux = zeros(E, length(dofb))
_interpolate_values_kernel!(
free_dofs,diri_dofs,uphys,celldofs,dofb,aux)
return free_dofs, diri_dofs




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
