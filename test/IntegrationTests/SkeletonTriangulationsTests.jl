module SkeletonTriangulationTests

using Gridap

model = CartesianDiscreteModel(partition=(3,4,2))

tags = [27,]
grid = SkeletonGrid(model,tags)

trian = Triangulation(grid)
descr1 = grid.descriptor1
descr2 = grid.descriptor2

strian = SkeletonTriangulation(trian,descr1,descr2)

test_triangulation(strian)

strian = SkeletonTriangulation(model)
#writevtk(strian,"strian")

end # module
