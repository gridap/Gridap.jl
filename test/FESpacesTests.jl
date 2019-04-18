##
using Numa, Test

using Numa.Quadratures
using Numa.Polytopes
using Numa.RefFEs
# using Numa.Meshes
using Numa.FESpaces
using Numa.FESpaces: ConformingFESpace

using Numa.Polytopes: PointInt
using Numa.CellValues

using Numa.Maps
using Numa.FieldValues

using Numa.Meshes

using Numa.CellValues: CellVectorByComposition

##
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
order=3
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
mesh = StructHexMesh(nparts)
fesp = ConformingFESpace{D,D,ScalarValue}(reffe,mesh)
using Numa.FESpaces: globaldofs
gldofs = globaldofs(fesp)
ndofs = gldofs[end][end]
@test gldofs[end][end]==(nparts1d*order+1)^D
cell_to_dofs = CellVectorByComposition(mesh.cellvefs, gldofs)
assembler = ConformingAssembler(fesp)
@test assembler.assembly_op_cols[1] == cell_to_dofs[1]
##
