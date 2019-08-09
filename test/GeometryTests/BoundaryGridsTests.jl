module BoundaryGridsTests

using Test
using Gridap

model = CartesianDiscreteModel(partition=(3,4,2))

tags = [23,24,25]
grid = BoundaryGrid(model,tags)
test_grid(grid,60,20)

tags = "boundary"
grid = BoundaryGrid(model,tags)

tags = ["physical_tag_23","physical_tag_24","physical_tag_25"]
grid = BoundaryGrid(model,tags)
test_grid(grid,60,20)

#writevtk(model,"model")

#writevtk(grid,"grid",celldata=["cell"=>grid.descriptor.facet_to_cell,"lfacet"=>grid.descriptor.facet_to_lfacet])

trian = Triangulation(grid)
@test isa(trian,BoundaryTriangulation)

end # module
