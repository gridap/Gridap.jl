module BoundaryTriangulationsTests

using Gridap

model = CartesianDiscreteModel(partition=(3,4,2))

tags = [23,24,25]
grid = BoundaryGrid(model,tags)

trian = Triangulation(grid)
descr = grid.descriptor

btrian = BoundaryTriangulation(trian,descr)

test_triangulation(btrian)

end # module
