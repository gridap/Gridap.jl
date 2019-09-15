module RaviartThomasFESpacesTests
##
using Test
using Gridap
using Gridap.MapApply
using Gridap.CellValuesGallery
using LinearAlgebra

D = 2
order = 1

dom = fill(2,D)
model = CartesianDiscreteModel(partition=tuple(dom...))

p = Polytope(fill(HEX_AXIS,D)...)
reffe = RaviartThomasRefFE(p,order)

grid = Grid(model,D)
trian = Triangulation(grid)
graph = GridGraph(model)
graph.data[2,1]

labels = FaceLabels(model)
# tags = [1,2,3,4]
tags = Int[]
fesp = ConformingFESpace(reffe,trian,graph,labels,tags)

fun(x) = VectorValue(1.0,0.0)

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

uh = interpolate(fesp,fun)
uh.free_dofs

a = [1;2;3;4]
b = Matrix(1.0I,4,4)
b[2,2] = -1.0
b = [1;-1;1;1]

ca = ConstantCellValue(a,4)
cb = ConstantCellValue(b,4)
cb*ca

# So, we must create a cell values with +1,-1
# @santiagobadia : Create the local_to_global_dofs and global_to_local_dofs

e = fun - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10
el2

##

















cellvec = integrate(inner(bh,uh),trian,quad)

bh = FEBasis(fespace)
uh = zero(fespace)


nfree = 48
ndiri = 0



cellmat = integrate(inner(bh,bh),trian,quad)

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

@test celldofids(fespace).data == 1:nfree




#ufun(x) = sin(x[1]) * cos(x[2])
#
#uh = interpolate(fespace,ufun)
#
#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh])


end # module
