module RTFESpacesTests
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
reffe = RTRefFE(p,Float64,order)

grid = Grid(model,D)
trian = Triangulation(grid)
graph = GridGraph(model)
graph.data[2,1]

labels = FaceLabels(model)
tags = [5,6,7,8]
fun(x) = VectorValue(x[1],x[2])
# tags = Int[]

fesp = ConformingFESpace(reffe,trian,graph,labels,tags)
dfesp = TrialFESpace(fesp,fun)


trian = Triangulation(model)
quad = CellQuadrature(trian,degree=2)

uh = interpolate(dfesp,fun)

e = fun - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

##

end # module
