
#
#  7 ---9--- 8 ---12-- 9
#
#  |         |         |
#  |         |         |
#  |         |         |
#  10   4    11   5    13
#  |         |         |
#  |         |         |
#  |         |         |
#
#  4 ---2--- 5 ---8--- 6
#
#  |         | \   3   |
#  |         |  \      |
#  |         |   \     |
#  3    1    4    6    7
#  |         |     \   |
#  |         |  2   \  |
#  |         |       \ |
#
#  1 ---1--- 2 ---5--- 3
#

                        #1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
cell_to_vertices_data = [1,2,4,5,2,3,5,3,6,5,4,5,7,8,5,6,8,9]
cell_to_vertices_ptrs = Int32[1,      5,    8,    11,     15,    19]
cell_to_vertices = Table(cell_to_vertices_data,cell_to_vertices_ptrs)

                       #1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
vertex_to_cells_data = [1,1,2,2,3,1,4,1,2,3,4,5,3,5,4,4,5,5]
vertex_to_cells_ptrs = Int32[1,2,  4,  6,  8,        13, 15,16, 18,19]
vertex_to_cells = Table(vertex_to_cells_data,vertex_to_cells_ptrs)

cell_to_faces_data = [1,2,3,4,5,4,6,7,6,8,2,9,10,11,8,12,11,13]
cell_to_faces_ptrs = Int32[1,5,8,11,15,19]
cell_to_faces = Table(cell_to_faces_data, cell_to_faces_ptrs)

quad_to_lface_to_lvertices = [[1,2],[3,4],[1,3],[2,4]]
quad_to_lface_to_ftype = Int8[1,1,1,1]

tri_to_lface_to_lvertices = [[1,2],[1,3],[2,3]]
tri_to_lface_to_ftype = Int8[2,2,2]

quad = Int8(1)
tri = Int8(2)
cell_to_ctype = [quad, tri, tri, quad, quad]

ctype_to_lface_to_lvertices = [
  quad_to_lface_to_lvertices, tri_to_lface_to_lvertices]

ctype_to_lface_to_ftype = [quad_to_lface_to_ftype, tri_to_lface_to_ftype]

face_to_ftype = Int8[1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1]

face_to_vertices_data = [1,2,4,5,1,4,2,5,2,3,3,5,3,6,6,5,7,8,4,7,5,8,8,9,6,9]
face_to_vertices_ptrs = Int32[1,3,5,7,9,11,13,15,17,19,21,23,25,27]
face_to_vertices = Table(face_to_vertices_data, face_to_vertices_ptrs)

face_to_isboundary = Bool[true, false, true, false, true, false, true, false, true, true, false, true, true]
vertex_to_isboundary = Bool[true, true, true, true, false, true, true, true, true]

