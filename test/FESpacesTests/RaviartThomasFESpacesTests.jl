module RaviartThomasFESpacesTests
##
using Test
using Gridap
using Gridap.MapApply
using Gridap.CellValuesGallery
using LinearAlgebra

D = 2
order = 1

dom = fill(4,D)
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

fun(x) = VectorValue(x[1],x[2])

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

uh = interpolate(fesp,fun)
uh.free_dofs

e = fun - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

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
