module BoundaryGridsTests

using Test
using Gridap

model = CartesianDiscreteModel(partition=(3,4,2))

tags = [23,24,25]
grid = BoundaryGrid(model,tags)
test_grid(grid,60,20)

#writevtk(grid,"grid",celldata=["cell"=>grid.facet_to_cell,"lfacet"=>grid.facet_to_lfacet])

end # module
