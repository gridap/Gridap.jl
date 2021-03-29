module PeriodicBCTests

using Gridap.Arrays
using Gridap.Geometry

using Test


## 2D mesh with periodic BC in x
model = CartesianDiscreteModel((0,1,0,1),(3,3);isperiodic=(true,false))

# Test cells to vertices
@test model.grid_topology.n_m_to_nface_to_mfaces[3,1] ==
  [[1, 2, 4, 5],
   [2, 3, 5, 6],
   [3, 1, 6, 4],
   [4, 5, 7, 8],
   [5, 6, 8, 9],
   [6, 4, 9, 7],
   [7, 8, 10, 11],
   [8, 9, 11, 12],
   [9, 7, 12, 10]]

# Test cells to faces
@test model.grid_topology.n_m_to_nface_to_mfaces[3,2] ==
  [[1, 2, 3, 4],
   [5, 6, 4, 7],
   [8, 9, 7, 3],
   [2, 10, 11, 12],
   [6, 13, 12, 14],
   [9, 15, 14, 11],
   [10, 16, 17, 18],
   [13, 19, 18, 20],
   [15, 21, 20, 17]]

# Test faces to vertices
@test model.grid_topology.n_m_to_nface_to_mfaces[2,1] ==
  [[1, 2], [4, 5], [1, 4], [2, 5], [2, 3], [5, 6], [3, 6], [3, 1], [6, 4],
   [7, 8], [4, 7], [5, 8], [8, 9], [6, 9], [9, 7],
   [10, 11], [7, 10], [8, 11], [11, 12], [9, 12], [12, 10]]

# Test vertices to faces
@test model.grid_topology.n_m_to_nface_to_mfaces[1,2] ==
  [[1, 3, 8],
   [1, 4, 5],
   [5, 7, 8],
   [2, 3, 9, 11],
   [2, 4, 6, 12],
   [6, 7, 9, 14],
   [10, 11, 15, 17],
   [10, 12, 13, 18],
   [13, 14, 15, 20],
   [16, 17, 21],
   [16, 18, 19],
   [19, 20, 21]]

# Test vertices to cells
@test model.grid_topology.n_m_to_nface_to_mfaces[1,3] ==
  [[1, 3],
   [1, 2],
   [2, 3],
   [1, 3, 4, 6],
   [1, 2, 4, 5],
   [2, 3, 5, 6],
   [4, 6, 7, 9],
   [4, 5, 7, 8],
   [5, 6, 8, 9],
   [7, 9],
   [7, 8],
   [8, 9]]

# Test faces to cells
@test model.grid_topology.n_m_to_nface_to_mfaces[2,3] ==
  [[1], [1, 4], [1, 3], [1, 2], [2], [2, 5], [2, 3], [3], [3, 6],
   [4, 7], [4, 6], [4, 5], [5, 8], [5, 6], [6, 9],
   [7], [7, 9], [7, 8], [8], [8, 9], [9]]

## 2D mesh with periodic BC in x & y
model = CartesianDiscreteModel((0,1,0,1),(3,3);isperiodic=(true,true))

 # Test cells to vertices
@test model.grid_topology.n_m_to_nface_to_mfaces[3,1] ==
  [[1, 2, 4, 5],
   [2, 3, 5, 6],
   [3, 1, 6, 4],
   [4, 5, 7, 8],
   [5, 6, 8, 9],
   [6, 4, 9, 7],
   [7, 8, 1, 2],
   [8, 9, 2, 3],
   [9, 7, 3, 1]]

# Test cells to faces
@test model.grid_topology.n_m_to_nface_to_mfaces[3,2] ==
  [[1, 2, 3, 4],
   [5, 6, 4, 7],
   [8, 9, 7, 3],
   [2, 10, 11, 12],
   [6, 13, 12, 14],
   [9, 15, 14, 11],
   [10, 1, 16, 17],
   [13, 5, 17, 18],
   [15, 8, 18, 16]]

# Test faces to vertices
@test model.grid_topology.n_m_to_nface_to_mfaces[2,1] ==
  [[1, 2], [4, 5], [1, 4], [2, 5], [2, 3], [5, 6], [3, 6], [3, 1], [6, 4],
   [7, 8], [4, 7], [5, 8], [8, 9], [6, 9], [9, 7],
   [7, 1], [8, 2], [9, 3]]

# Test vertices to faces
@test model.grid_topology.n_m_to_nface_to_mfaces[1,2] ==
  [[1, 3, 8, 16],
   [1, 4, 5, 17],
   [5, 7, 8, 18],
   [2, 3, 9, 11],
   [2, 4, 6, 12],
   [6, 7, 9, 14],
   [10, 11, 15, 16],
   [10, 12, 13, 17],
   [13, 14, 15, 18]]

# Test vertices to cells
@test model.grid_topology.n_m_to_nface_to_mfaces[1,3] ==
  [[1, 3, 7, 9],
   [1, 2, 7, 8],
   [2, 3, 8, 9],
   [1, 3, 4, 6],
   [1, 2, 4, 5],
   [2, 3, 5, 6],
   [4, 6, 7, 9],
   [4, 5, 7, 8],
   [5, 6, 8, 9]]

# Test faces to cells
@test model.grid_topology.n_m_to_nface_to_mfaces[2,3] ==
  [[1, 7], [1, 4], [1, 3], [1, 2], [2, 8], [2, 5], [2, 3], [3, 9], [3, 6],
   [4, 7], [4, 6], [4, 5], [5, 8], [5, 6], [6, 9],
   [7, 9], [7, 8], [8, 9]]

end #module
