##
using Test

using Gridap
using Gridap.Quadratures
using Gridap.CellQuadratures
# using Gridap.CellIntegration
# using Gridap.CellValues
using Gridap.CellMaps
import Gridap.CellIntegration: cellcoordinates, cellbasis

using Gridap.CellValues.ConstantCellValues
using Gridap.CellValues: IndexCellArray
import Gridap.CellValues: cellsize
using Gridap.Polytopes
using Gridap.Polytopes: PointInt
using Gridap.RefFEs
using Gridap.FieldValues

import Gridap: gradient

using Gridap.Meshes
using Gridap.FESpaces: ConformingFESpace

using Gridap.CellIntegration

using Gridap.FESpaces: ConformingAssembler
using SparseArrays
using Gridap.FESpaces: Assembler

using Gridap.Geometry
using Gridap.Geometry.Cartesian

##
D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = Triangulation(grid) # Generates the Triangulation associated with this grid
graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
phi = geomap(trian)
typeof(phi)
order=1
orders=order*ones(Int64,D)
pol_array = celltypes(trian)
extrusion = PointInt{D}(pol_array[1])
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = ConstantCellValue(basis, ncells(trian))
quad = quadrature(trian,order=2)

fun(x::Point{2}) = x[1]^2
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
uphys = fun ∘ phi
typeof(uphys)
qps = coordinates(quad)
evaluate(uphys,qps)
nods = [Point{2}(-1.0,-1.0), Point{2}(1.0,-1.0), Point{2}(-1.0,1.0), Point{2}(1.0,1.0)]
cnods = ConstantCellValue(nods,ncells(trian))
evaluate(phi,cnods)
evaluate(uphys,cnods)
##

using Gridap.FESpaces
using Gridap.FESpaces: interpolate
fesp = ConformingFESpace(reffe,trian,graph)
assembler = ConformingAssembler(fesp)
funh = interpolate(fun, fesp)

# @santiagobadia : To do, replace RefFE with an array of RefFEs!!!

##
v = zeros(Int64, assembler.num_dofs)
for l2g in celldofs
  v[l2g] = l2g
end
@show v
@test v == [ i for i in 1:assembler.num_dofs]


# Zero: Here I want to extract the gdofs of the vefs of a cell and create a global
# vector using the Assembler, of the global dofs. The test would be gid[i] = i
# In a second step, do it for Dirichlet too.

# Two: Compute the free and dof vectors for an analytical function, see below

# Three: Create the FEFunction







Interpolator(fesp, fun, phi)
# @santiagobadia : Don't we want a triangulation not a grid in fespace



isa(uphys, CellField)
isa(dofb, CellValue)
using Gridap.CellValues
struct Interpolator{M<:CellField,V<:DOFBasis}
  u::M
  dofb::V
  # cached vector, etc...
  # Probably better just CellValue{<:RefFE} from FESpace
  # How to make it efficient ?
end
using Gridap.FESpaces
function Interpolator(this::FESpace, fun::Function, ass::Assembler)
  reffe = this.reffe
  dofb = reffe.dofbasis
  trian = fesp.trian
  phi = geomap(trian)
  uphys = fun ∘ phi
  celldofs = assembler.assembly_op_rows
  v = zeros(Int64, assembler.num_dofs)
  for (imap,l2g) in zip(uphys,celldofs)
      dofs[l2g] = evaluate(dofb,imap)
    end
end

  # Interpolator(uphys, dofb)
end

function iterate
end

function computevals!
  evaluate(dofb, cmap)
end

##




basis = cellbasis(trian)
physbasis = attachgeomap(basis,phi)
points = (trian)
ab(v,u) = inner(∇(v),∇(u)) #+ inner(v,u)
V = physbasis
U = physbasis
# fun(x::Point{2}) = x[1]*x[2] + x[1]
# gradfun(x::Point{2}) = VectorValue(x[2] + 1.0, x[1])
ksca = integrate(ab(uphys,uphys),trian,quad)
sum(ksca)
kvec = integrate(ab(V,uphys),trian,quad)
typeof(kvec)
kmat = integrate(ab(V,U),trian,quad)
##


for (i,cell) in enumerate(trian)
  @show cell
  @show i
end




# @santiagobadia : Not tested yet
using Gridap.FESpaces: assemble
sys_vec = assemble(assembler,kvec)
sys_mat = assemble(assembler,kmat)
