module GridTopologiesTests

using Test
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: GridTopologyMock

function compute_cell_perm_indices(top::GridTopology,d::Integer)

  D = num_cell_dims(top)
  face_to_fvertex_to_vertex = Table(get_faces(top,d,0))
  face_to_ftype = get_face_reffe_type(top,d)
  reffaces = get_reffes(NodalReferenceFE{d},top)
  facepolytopes = map(get_polytope,reffaces)
  ftype_to_pindex_to_cfvertex_to_fvertex = map(get_vertex_permutations,facepolytopes)
  cell_to_cvertex_to_vertex = Table(get_faces(top,D,0))
  cell_to_lface_to_face = Table(get_faces(top,D,d))
  cell_to_ctype = get_cell_type(top)
  reffes = get_reffes(top)
  polytopes = map(get_polytope,reffes)
  ctype_to_lface_to_cvertices = map( (p)->get_faces(p,d,0), polytopes )
  data = similar(cell_to_lface_to_face.data)
  ptrs = cell_to_lface_to_face.ptrs
  cell_to_lface_to_pindex = Table(data,ptrs)

  _compute_cell_perm_indices!(
    cell_to_lface_to_pindex,
    cell_to_lface_to_face,
    cell_to_cvertex_to_vertex,
    cell_to_ctype,
    ctype_to_lface_to_cvertices,
    face_to_fvertex_to_vertex,
    face_to_ftype,
    ftype_to_pindex_to_cfvertex_to_fvertex)

  cell_to_lface_to_pindex
end

top = GridTopologyMock()
test_grid_topology(top)

@test num_dims(top) == 2
@test num_cell_dims(top) == 2
@test num_point_dims(top) == 2
@test num_faces(top) == 9 + 13 + 5
@test num_faces(top,0) == 9
@test num_faces(top,1) == 13
@test num_faces(top,2) == 5
@test num_cells(top) == 5
@test num_vertices(top) == 9
@test num_edges(top) == 13
@test num_facets(top) == 13
@test get_dimranges(top) == [1:9, 10:22, 23:27]
@test get_dimrange(top,0) == 1:9
@test get_dimrange(top,1) == 10:22
@test get_dimrange(top,2) == 23:27
@test get_offsets(top) == [0, 9, 22]
@test get_offset(top,0) == 0
@test get_offset(top,1) == 9
@test get_offset(top,2) == 22

cell_to_faces = get_cell_faces(top)
@test cell_to_faces == [
  [1,2,4,5,10,11,12,13,23],[2,3,5,14,13,15,24],[3,6,5,16,15,17,25],
  [4,5,7,8,11,18,19,20,26],[5,6,8,9,17,21,20,22,27]]
@test compute_cell_faces(top) == get_cell_faces(top)

@test get_face_vertices(top,0) == get_faces(top,0,0)
@test get_face_vertices(top,1) == get_faces(top,1,0)
@test get_face_vertices(top,2) == get_faces(top,2,0)

face_to_vertices = get_face_vertices(top)
@test face_to_vertices == vcat( [ get_face_vertices(top,d) for d in 0:2]... )
@test compute_face_vertices(top) == get_face_vertices(top)

@test is_simplex(top) == false
@test is_n_cube(top) == false

reffaces, face_to_ftype = compute_reffaces(Polytope{1},top)
@test reffaces == [SEGMENT,]
@test face_to_ftype == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
@test reffaces == get_reffaces(Polytope{1},top)
@test face_to_ftype == get_face_type(top,1)

@test compute_isboundary_face(top, 0) == Bool[1, 1, 1, 1, 0, 1, 1, 1, 1]
@test compute_isboundary_face(top, 1) == Bool[1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1]
@test compute_isboundary_face(top, 2) == Bool[1, 1, 1, 1, 1]

@test get_isboundary_face(top, 0) == compute_isboundary_face(top, 0)
@test get_isboundary_face(top, 1) == compute_isboundary_face(top, 1)
@test get_isboundary_face(top, 2) == compute_isboundary_face(top, 2)

@test compute_isboundary_face(top) == Bool[1,1,1,1,0,1,1,1,1,1,0,1,0,1,0,1,0,1,1,0,1,1,1,1,1,1,1]
@test compute_isboundary_face(top) == get_isboundary_face(top)

end # module
