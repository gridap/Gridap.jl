module SkeletonGridsTests

using Test
using Gridap

model = CartesianDiscreteModel(partition=(5,4,3))

tags = [27,]
grid = SkeletonGrid(model,tags)
test_grid(grid,120,133)

#writevtk(grid,"grid",celldata=[
#  "cell1"=>grid.descriptor1.facet_to_cell,
#  "lfacet1"=>grid.descriptor1.facet_to_lfacet,
#  "cell2"=>grid.descriptor2.facet_to_cell,
#  "lfacet2"=>grid.descriptor2.facet_to_lfacet ])

trian = Triangulation(grid)
@test isa(trian,SkeletonTriangulation)

end # module
