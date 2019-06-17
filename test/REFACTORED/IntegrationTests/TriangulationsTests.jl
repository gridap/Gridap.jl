module TriangulationsTests

using Gridap
using Test

partition = (2,4)
trian = CartesianTriangulation(partition)

test_triangulation(trian)

@test num_cells(trian) == 8

coords = CellPoints(trian)

x4 = Point{2,Float64}[[0.0, -0.5], [1.0, -0.5], [0.0, 0.0], [1.0, 0.0]] 

@test coords[4] ≈ x4

quad = CellQuadrature(trian,order=2)

end # module
