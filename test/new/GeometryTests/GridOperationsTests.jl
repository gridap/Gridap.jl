module GridOperationsTests

using Test

using Gridap.Helpers
using Gridap.Arrays
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

end # module
