module GridOperationsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry


using Gridap.Helpers
using Gridap.Arrays
using Gridap.ReferenceFEs: _find_unique_with_indices
include("../../../src/new/Geometry/GridOperations.jl")

include("Mock2d.jl")

_vertex_to_cells = generate_cells_around(cell_to_vertices)
@test isa(_vertex_to_cells,Table{Int,Int32})
@test _vertex_to_cells == vertex_to_cells

_cell_to_faces = generate_cell_to_faces(
  cell_to_vertices,
  ctype_to_lface_to_lvertices,
  cell_to_ctype,
  vertex_to_cells)

@test isa(_cell_to_faces,Table{Int,Int32})
@test _cell_to_faces == cell_to_faces

_face_to_ftype = generate_face_to_face_type(
  cell_to_faces,
  cell_to_ctype,
  ctype_to_lface_to_ftype)

@test isa(_face_to_ftype,Vector{Int8})
@test face_to_ftype == _face_to_ftype

_face_to_vertices = generate_face_to_vertices(
  cell_to_vertices,
  cell_to_faces,
  cell_to_ctype,
  ctype_to_lface_to_lvertices)

@test isa(_face_to_vertices,Table{Int,Int32})
@test _face_to_vertices == face_to_vertices

vertex_to_faces = generate_cells_around(face_to_vertices)

_cell_to_faces = find_cell_to_faces(
  cell_to_vertices,
  ctype_to_lface_to_lvertices,
  cell_to_ctype,
  vertex_to_faces)

@test isa(_cell_to_faces,Table{Int,Int32})
@test _cell_to_faces == cell_to_faces

face_to_cells = generate_cells_around(cell_to_faces)

_face_to_isboundary = generate_facet_to_isboundary(face_to_cells)
@test _face_to_isboundary == face_to_isboundary

_vertex_to_isboundary = generate_face_to_isboundary(face_to_isboundary,vertex_to_faces)
@test _vertex_to_isboundary == _vertex_to_isboundary

cell_to_nodes = get_cell_nodes(grid)
cell_to_type = get_cell_types(grid)
type_to_reffes = get_reffes(grid)

type_to_vertex_to_node = map(get_vertex_node, type_to_reffes)

_cell_to_vertices, vertex_to_node = _generate_cell_to_vertices(
  cell_to_nodes, cell_to_ctype, type_to_vertex_to_node)

@test isa(_cell_to_vertices,Table{Int,Int32})

_cell_to_faces = generate_cell_to_faces(1,grid,cell_to_vertices,vertex_to_cells)
@test _cell_to_faces == cell_to_faces
@test isa(_cell_to_faces,Table{Int,Int32})


@test vertex_to_node == [1, 2, 4, 5, 3, 6, 7, 8, 9]
@test _cell_to_vertices == [[1, 2, 3, 4], [2, 5, 4], [5, 6, 4], [3, 4, 7, 8], [4, 6, 8, 9]]

_cell_to_vertices, vertex_to_node = generate_cell_to_vertices(grid)

@test _cell_to_vertices == cell_to_vertices
@test isa(_cell_to_vertices,Table{Int,Int32})

ctype_to_reffe = get_reffes(grid)
ftype_to_refface, ctype_to_lface_to_ftype = _generate_ftype_to_refface(Val{1}(),ctype_to_reffe)

@test length(ftype_to_refface) == 1
@test ctype_to_lface_to_ftype == [[1, 1, 1, 1], [1, 1, 1]]

face_to_ftype =  generate_face_to_face_type(cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype)
@test length(face_to_ftype) == length(face_to_vertices)
@test face_to_ftype == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

end # module
