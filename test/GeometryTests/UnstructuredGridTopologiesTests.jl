module UnstructuredGridTopologiesTests

using Test
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: GridTopologyMock

m = GridTopologyMock()

g = UnstructuredGridTopology(
  get_vertex_coordinates(m),
  get_cell_vertices(m),
  get_cell_type(m),
  get_polytopes(m))

test_grid_topology(g)
test_grid_topology(g)

@test num_faces(g,0) == num_faces(m,0)
@test num_faces(g,1) == num_faces(m,1)
@test num_faces(g,2) == num_faces(m,2)
@test get_reffaces(Polytope{1},g) == get_reffaces(Polytope{1},m)
@test get_reffaces(Polytope{0},g) == get_reffaces(Polytope{0},m)
@test get_reffaces(Polytope{2},g) == get_reffaces(Polytope{2},m)
@test get_face_type(g,1) == get_face_type(m,1)
@test get_isboundary_face(g,0) == get_isboundary_face(m,0)
@test get_isboundary_face(g,1) == get_isboundary_face(m,1)
@test get_isboundary_face(g,2) == get_isboundary_face(m,2)
@test get_faces(g,2,0) == get_faces(m,2,0)
@test get_faces(g,2,1) == get_faces(m,2,1)
@test get_faces(g,1,2) == get_faces(m,1,2)
@test get_faces(g,2,1) == get_faces(m,2,1)
@test get_faces(g,0,1) == get_faces(m,0,1)
@test get_faces(g,1,0) == get_faces(m,1,0)
@test get_faces(g,1,1) == get_faces(m,1,1)
@test get_vertex_coordinates(g) == get_vertex_coordinates(m)           
@test is_oriented(g) == false

g = UnstructuredGridTopology(
  get_vertex_coordinates(m),
  get_cell_vertices(m),
  get_cell_type(m),
  get_polytopes(m),
  Oriented())

@test is_oriented(g) == true

g = UnstructuredGridTopology(
  get_vertex_coordinates(m),
  [get_face_vertices(m,d) for d in 0:num_cell_dims(m)],
  get_cell_type(m),
  get_polytopes(m))

@test get_faces(g,2,0) == get_faces(m,2,0)
@test get_faces(g,2,1) == get_faces(m,2,1)
@test get_faces(g,1,2) == get_faces(m,1,2)
@test get_faces(g,2,1) == get_faces(m,2,1)
@test get_faces(g,0,1) == get_faces(m,0,1)
@test get_faces(g,1,0) == get_faces(m,1,0)
@test get_faces(g,1,1) == get_faces(m,1,1)
@test is_oriented(g) == false

test_grid_topology(g)
test_grid_topology(g)

g = UnstructuredGridTopology(
  get_vertex_coordinates(m),
  [get_face_vertices(m,d) for d in 0:num_cell_dims(m)],
  get_cell_type(m),
  get_polytopes(m),
  Oriented())

@test is_oriented(g) == true

m = GridTopologyMock()
topo = UnstructuredGridTopology(m)
test_grid_topology(topo)
@test topo === UnstructuredGridTopology(topo)

# Extract grid topology

grid = GridMock()
topo = GridTopology(grid)
test_grid_topology(topo)

# Extract grid topology

domain = (1,2,1,2)
partition = (3,3)
grid = CartesianGrid(domain,partition)
topo = GridTopology(grid)
test_grid_topology(topo)

end # module
