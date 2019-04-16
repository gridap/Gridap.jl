##
using Numa
using Numa.Quadratures
using Numa.CellQuadratures
# using Numa.CellIntegration
# using Numa.CellValues
using Numa.CellFunctions
import Numa.CellIntegration: cellcoordinates, cellbasis

using Numa.CellValues: IndexCellArray
import Numa.CellValues: cellsize
using Numa.Polytopes
using Numa.Polytopes: PointInt
using Numa.RefFEs
using Numa.FieldValues

import Numa: gradient
using Numa.CellValues: ConstantCellValue

using Numa.Meshes
using Numa.FESpaces: ConformingFESpace

using Numa.CellIntegration

using Numa.FESpaces: ConformingAssembler
using SparseArrays
using Numa.FESpaces: Assembler

##
D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
phi = geomap(trian)
order=1
orders=order*ones(Int64,D)
pol_array = celltypes(trian)
pol_array[1]
D
extrusion = PointInt{D}(pol_array[1])
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = CellBasisFromSingleInterpolation(basis)
quad = quadrature(trian,order=2)
phi = geomap(trian)
basis = cellbasis(trian)
physbasis = attachgeomap(basis,phi)
ab(v,u) = inner(∇(v),∇(u)) #+ inner(v,u)
V = physbasis
U = physbasis
# fun(x::Point{2}) = x[1]*x[2] + x[1]
# gradfun(x::Point{2}) = VectorValue(x[2] + 1.0, x[1])
fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
uphys = fun ∘ phi
ksca = integrate(ab(uphys,uphys),trian,quad)
sum(ksca)
kvec = integrate(ab(V,uphys),trian,quad)
typeof(kvec)
kmat = integrate(ab(V,U),trian,quad)
##
fesp = ConformingFESpace{D,D,ScalarValue}(reffe,grid)
assembler = ConformingAssembler(fesp)
# @santiagobadia : Not tested yet
using Numa.FESpaces: assemble
sys_vec = assemble(assembler,kvec)
sys_mat = assemble(assembler,kmat)
