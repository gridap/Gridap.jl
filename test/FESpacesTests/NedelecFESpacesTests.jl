module NedelecFESpacesTests
##
using Test
using Gridap
using Gridap.MapApply
using Gridap.CellValuesGallery
using LinearAlgebra

D = 2
order = 1

dom = fill(1,D)
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(1,1))

p = Polytope(fill(HEX_AXIS,D)...)
reffe = NedelecRefFE(p,order)

grid = Grid(model,D)
trian = Triangulation(grid)
graph = GridGraph(model)



labels = FaceLabels(model)
# tags = [5,6,7,8]
tags = Int[]
fun(x) = VectorValue(x[1],x[2])
# fun(x) = VectorValue(1.0,1.0)
# tags = Int[]

fesp = ConformingFESpace(reffe,trian,graph,labels,tags)
dfesp = TrialFESpace(fesp,fun)


trian = Triangulation(model)
quad = CellQuadrature(trian,degree=2)

uh = interpolate(dfesp,fun)
uh.free_dofs
uh.diri_dofs

uh2 = FEFunction(uh.fespace,0.5*uh.free_dofs,uh.diri_dofs)

e = fun - uh
uh.free_dofs
uh.diri_dofs

ps = ConstantCellVector(1,1)
isa(ps,CellValue)
evaluate(uh,ps)
length(uh)
length(ps)



el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

##

end # module
