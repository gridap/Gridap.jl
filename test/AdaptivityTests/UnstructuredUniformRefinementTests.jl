module UnstructuredUniformRefinementTests
  
using Test
using Gridap
using Gridap.Adaptivity
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Adaptivity: test_unstructured_uniform_refinement
using Gridap.Adaptivity: cube_simplex_pattern_dimfid
using Gridap.Adaptivity: cube_simplex_reference_grid
using Gridap.Adaptivity: compute_d_dface_offsets
using Gridap.Adaptivity: compute_cell_offsets
using Gridap.Adaptivity: unstructured_uniform_cell_l2gmap_and_nnodes
using Gridap.Adaptivity: unstructured_uniform_coordinates
using Gridap.Adaptivity: unstructured_uniform_connectivity
using Gridap.Adaptivity: unstructured_uniform_topology

test_unstructured_uniform_refinement()

#### TRI ####

n = 2
p = TRI

face_bary = VectorValue{2,Float64}[
  (0,0),(1,0),(0,1),(0.5,0),(0,0.5),(0.5,0.5),(1/3,1/3)
]
@test face_bary == map(mean,get_face_coordinates(p))

# 0->-1, 1->0, o->1
# sum_one->1
pt_dimfid = Dict(
  # coordinates (0,0), dimfid (0,1), pattern (-1,-1,0)
  (-1,-1,0) => (0,1),
  # coordinates (0.5,0.5), dimfid (1,3), pattern (1,1,1)
  (1,1,1) => (1,3),
  # (1,0), (0,2), (0,-1,1)
  (0,-1,1) => (0,2),
  # (0,1), (0,3), (-1,0,1)
  (-1,0,1) => (0,3),
  # (0.5,0), (1,1), (1,-1,0)
  (1,-1,0) => (1,1),
  # (0,0.5), (1,2), (-1,1,0)
  (-1,1,0) => (1,2),
  # (1/3,1/3), (2,1), (1,1,0)
  (1,1,0) => (2,1)
)
pt_map = cube_simplex_pattern_dimfid(p)
@test all( [pt_map[k] == pt_dimfid[k] for k in keys(pt_map)] )

# 3
# | \
# |  \
# 5 - 6
# | \ |\
# |  \| \
# 1 - 4 - 2
ref_grid = cube_simplex_reference_grid(p,n,pt_map)
coords = VectorValue{2,Float64}[
  (0,0),(1,0),(0,1),(0.5,0),(0.,0.5),(0.5,0.5)
]
cell_ids = [
  [1,4,5], [5,6,3], [4,2,6], [4,6,5]
]

test_grid(ref_grid)
@test coords == get_node_coordinates(ref_grid)
@test cell_ids == get_cell_node_ids(ref_grid)


#### QUAD ####

p = QUAD
n = 2

face_bary = VectorValue{2,Float64}[
  (0,0),(1,0),(0,1),(1,1),
  (0.5,0),(0.5,1),(0,0.5),(1,0.5),
  (0.5,0.5)
]

@test face_bary == map(mean,get_face_coordinates(p))

# 0->-1, 1->0, o->1
pt_dimfid = Dict(
  (-1,-1) => (0,1),
  (0,-1) => (0,2),
  (-1,0) => (0,3),
  (0,0) => (0,4),
  (1,-1) => (1,1),
  (1,0) => (1,2),
  (-1,1) => (1,3),
  (0,1) => (1,4),
  (1,1) => (2,1)
)
pt_map = cube_simplex_pattern_dimfid(p)
@test all( [pt_map[k] == pt_dimfid[k] for k in keys(pt_map)] )


# 3 ---- 6 ---- 4
# |      |      |
# |      |      |
# 7 ---- 9 ---- 8
# |      |      |
# |      |      |
# 1 ---- 5 ---- 2
ref_grid = cube_simplex_reference_grid(p,n,pt_map)
coords = VectorValue{2,Float64}[
  (0,0),(1,0),(0,1),(1,1),
  (0.5,0),(0.5,1),(0,0.5),(1,0.5),
  (0.5,0.5)
]
cell_ids = [
  [1,5,7,9], [5,2,9,8],
  [7,9,3,6], [9,8,6,4]
]

test_grid(ref_grid)
@test coords == get_node_coordinates(ref_grid)
@test cell_ids == get_cell_node_ids(ref_grid)

end