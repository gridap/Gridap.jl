module GridGraphTests

using Test

using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: GridGraphMock

g = GridGraphMock()
test_grid_graph(g)

@test num_dims(g) == 2
@test num_cells(g) == 5
@test num_faces(g) == 27
@test num_facets(g) == 13
@test num_edges(g) == 13
@test num_vertices(g) == 9
@test get_facedims(g) == [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2]
@test get_offsets(g) == [0,9,22]
@test get_offset(g,0) == 0
@test get_offset(g,1) == 9
@test get_offset(g,2) == 22

end # module
