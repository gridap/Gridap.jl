
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

struct ConformingTrianMock <: ConformingTriangulation{2,2} end

function get_node_coordinates(::ConformingTrianMock)
  Point{2,Float64}[(0,0),(1,0),(2,0),(1,1),(2,1),(0,2),(2,2)]
end

function get_cell_nodes(::ConformingTrianMock)
                          #1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
  cell_to_vertices_data = [1,2,4,5,2,3,5,3,6,5,4,5,7,8,5,6,8,9]
  cell_to_vertices_ptrs = Int32[1,      5,    8,    11,     15,    19]
  Table(cell_to_vertices_data,cell_to_vertices_ptrs)
end

function get_reffes(::ConformingTrianMock)
  [QUAD4, TRI3]
end

function get_cell_types(::ConformingTrianMock)
  quad = Int8(1)
  tri = Int8(2)
  [quad, tri, tri, quad, quad]
end

struct GridGraphMock <: GridGraph end

function num_dims(g::GridGraphMock)
  2
end

function get_dimranges(g::GridGraphMock)
  [1:9, 10:22, 23:27]
end

function get_cell_faces(g::GridGraphMock)
  t = Table([
   [1, 2, 4, 5, 10, 11, 12, 13, 23],
   [2, 3, 5, 14, 13, 15, 24],
   [3, 6, 5, 16, 15, 17, 25],
   [4, 5, 7, 8, 11, 18, 19, 20, 26],
   [5, 6, 8, 9, 17, 21, 20, 22, 27] ])
  convert(Table{Int,Int32},t)
end

function get_face_cells(g::GridGraphMock)
  face_to_cells = collect(get_face_cells(g,0))
  append!(face_to_cells,collect(get_face_cells(g,1)))
  append!(face_to_cells,collect(get_face_cells(g,2)))
  Table(face_to_cells)
end

function get_cell_faces(g::GridGraphMock,facedim::Integer)
  if facedim == 0
    cell_to_vertices_data = [1,2,4,5,2,3,5,3,6,5,4,5,7,8,5,6,8,9]
    cell_to_vertices_ptrs = Int32[1,5,8,11,15,19]
    return Table(cell_to_vertices_data,cell_to_vertices_ptrs)
  elseif facedim == 1
    cell_to_faces_data = [1,2,3,4,5,4,6,7,6,8,2,9,10,11,8,12,11,13]
    cell_to_faces_ptrs = Int32[1,5,8,11,15,19]
    return Table(cell_to_faces_data, cell_to_faces_ptrs)
  elseif facedim == 2
    ncells = 5
    return identity_table(Int,Int32,ncells)
  else
    @unreachable
  end
end

function get_face_cells(g::GridGraphMock,facedim::Integer)
  if facedim == 0
    vertex_to_cells_data = [1,1,2,2,3,1,4,1,2,3,4,5,3,5,4,4,5,5]
    vertex_to_cells_ptrs = Int32[1,2,4,6,8,13,15,16,18,19]
    return Table(vertex_to_cells_data,vertex_to_cells_ptrs)
  elseif facedim == 1
    face_to_cells_data = [1,1,4,1,1,2,2,2,3,3,3,5,4,4,4,5,5,5]
    face_to_cells_ptrs = Int32[1,2,4,5,7,8,10,11,13,14,15,17,18,19]
    return Table(face_to_cells_data,face_to_cells_ptrs)
  elseif facedim == 2
    ncells = 5
    return identity_table(Int,Int32,ncells)
  else
    @unreachable
  end
end

function get_vertex_node(g::GridGraphMock)
  collect(1:9)
end

