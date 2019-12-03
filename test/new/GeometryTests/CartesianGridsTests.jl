module CartesianGridsTests

using FillArrays
using Test

using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Geometry
using Gridap.ReferenceFEs

domain = (0.0,1.0,-1.0,2.0)
partition = (3,4)

grid = CartesianGrid(domain,partition)
test_conforming_triangulation(grid)

@test num_cell_dims(grid) == 2
@test num_point_dims(grid) == 2
@test num_cells(grid) == 3*4
@test num_nodes(grid) == (3+1)*(4+1)

x = get_node_coordinates(grid)
@test length(x) == 20
@test x[13] == Point(0.0,1.25)
@test x[4] == Point(1.0,-1.0)
@test x[1] == Point(0.0,-1.0)
@test x[end] == Point(1.0,2.0)

t = get_cell_nodes(grid)
@test t[1] == [1,2,5,6]
@test t[11] == [14,15,18,19]

cell_type = get_cell_type(grid)
@test isa(cell_type,Fill)

reffes = get_reffes(grid)
@test length(reffes) == 1
@test reffes[1] == QUAD4

ugrid = UnstructuredGrid(grid)
test_conforming_triangulation(ugrid)

end # module
