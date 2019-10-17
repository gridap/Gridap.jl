module RefFEsTests
##
using Gridap.DOFBases: LagrangianDOFBasis

using Gridap, Test

using Gridap.Helpers

using Base.Cartesian

import Gridap: ∇
# 1D reffe

F = Gridap.ConformingFESpaces

D = 1
model = CartesianDiscreteModel(partition=(4,))
diri = "boundary"
p = Polytope(1,)
order = 4

function MultiShapeConformingFESpace(
  reffe::RefFE{D,T},
  trian::Triangulation{D,Z},
  graph::GridGraph,
  labels::FaceLabels,
  diri_tags::Vector{Int}) where {D,Z,T}
  args = _setup_multishape_conforming_fe_fields(reffe,trian,graph,labels,diri_tags,D)
  ConformingFESpace{D,Z,T}(args...)
end

function _setup_multishape_conforming_fe_fields(reffe,trian,graph,labels,diri_tags,D)
  dim_to_nface_eqclass, nfree, ndiri  = F._generate_dim_to_nface_to_dofs(
    reffe, graph, labels, diri_tags)
  cellvefs_dim = [connections(graph,D,i) for i in 0:D]
  offset = length.(dim_to_nface_eqclass)
  for i in 2:length(offset)
    offset[i] += offset[i-1]
  end
  offset = tuple(offset[1:(end-1)]...)
  cellvefs = local_append(offset, cellvefs_dim...)
  dofs_all = append(dim_to_nface_eqclass...)
  cell_eqclass = F.CellEqClass(cellvefs, dofs_all, reffe)
  values = [shfbasis(reffe)]
  pointer = [ 1 for i in 1:ncells(trian)]
  shb = CompressedCellValue(values,pointer)
  shb = ConstantCellValue(shfbasis(reffe), ncells(trian))
  phi = CellGeomap(trian)
  basis = attachgeomap(shb,phi)
  return dim_to_nface_eqclass, cell_eqclass, nfree, ndiri, diri_tags,
    reffe, trian, graph, labels, basis
end


reffe = LagrangianRefFE(Float64,p,order)
trian = Triangulation(model)
grid = Grid(model,D)
trian = Triangulation(grid)
orders = fill(order,D)
graph = GridGraph(model)
labels = FaceLabels(model)
mydiri_tags = Gridap.BoundaryGrids._setup_tags(labels,diri)
fespace = ConformingFESpace(reffe,trian,graph,labels,mydiri_tags)
# fespace = MultiShapeConformingFESpace(reffe,trian,graph,labels,mydiri_tags)


##

# Define manufactured functions
ufun(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = 0.0

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

end # module
