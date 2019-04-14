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

using Numa.Fields
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
mesh.cellvefs
cell_to_dofs = CellVectorByComposition(mesh.cellvefs, gldofs)
assembler = ConformingAssembler(fesp)
@test assembler.assembly_op_cols[1] == cell_to_dofs[1]
##








##
# Vector FE space
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
order=2
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE(polytope,orders)
gps=[2,2]
quad=Quadratures.TensorProductQuadrature(orders = tuple(2*orders...))
shfscal = shfsps(reffe,quad.coords)
reffe = LagrangianRefFE(polytope,orders,2)
shftens = shfsps(reffe,quad.coords)
@test shftens[36,1,2,2]==shfscal[9,1,1]
##
points=quad.coords
reffe = LagrangianRefFE(polytope,orders,0)
gradscal = gradshfsps(reffe,points)
reffe = LagrangianRefFE(polytope,orders,1)
gradtens = gradshfsps(reffe,points)
@test gradtens[10,1,2,2]==gradscal[1,1,2]
@test gradtens[10,3,2,2]==gradscal[1,3,2]
@test gradtens[1,3,2,1]==gradscal[1,3,2]
@test gradtens[1,3,2,2]==0.0
@test gradtens[10,3,2,1]==0.0
##
