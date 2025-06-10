
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

"""
    struct GridMock <: Grid{2,2}

For tests.
"""
struct GridMock <: Grid{2,2} end

function get_node_coordinates(::GridMock)
  Point{2,Float64}[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
end

function get_cell_node_ids(::GridMock)
                          #1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
  cell_to_vertices_data = [1,2,4,5,2,3,5,3,6,5,4,5,7,8,5,6,8,9]
  cell_to_vertices_ptrs = Int32[1,      5,    8,    11,     15,    19]
  Table(cell_to_vertices_data,cell_to_vertices_ptrs)
end

function get_reffes(::GridMock)
  [QUAD4, TRI3]
end

function get_cell_type(::GridMock)
  quad = Int8(1)
  tri = Int8(2)
  [quad, tri, tri, quad, quad]
end
