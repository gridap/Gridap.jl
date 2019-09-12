module BoundaryTriangulationsTests

using Gridap

model = CartesianDiscreteModel(partition=(3,4,2))

tags = [23,24,25]
grid = BoundaryGrid(model,tags)

trian = Triangulation(grid)
descr = grid.descriptor

btrian = BoundaryTriangulation(trian,descr)
test_triangulation(btrian)

tags = [23,24,25]
btrian = BoundaryTriangulation(model,tags)

tags = "boundary"
btrian = BoundaryTriangulation(model,tags)

tags = ["physical_tag_23","physical_tag_24","physical_tag_25"]
btrian = BoundaryTriangulation(model,tags)

btrian = BoundaryTriangulation(model)


end # module
