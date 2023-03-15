module RefinementRuleBoundaryTests

using Test
using Gridap
using Gridap.Adaptivity
using Gridap.Geometry
using Gridap.ReferenceFEs

function _get_terms(poly::Polytope,orders)
  _nodes, facenodes = Gridap.ReferenceFEs._compute_nodes(poly,orders)
  terms = Gridap.ReferenceFEs._coords_to_terms(_nodes,orders)
  return terms
end

poly = QUAD
rr = RefinementRule(QUAD,(2,2))

grid = rr.ref_grid
trian = Triangulation(grid)
btrian = Boundary(grid)

# This is what should be implemented, probably hardcoded for each refinement rule...
# For P4est, we would only need this one and the one for 3D 
c2f_edges = [[1,5],[8,11],[3,9],[7,12]]

coarse_reffe = ReferenceFE(QUAD,lagrangian,Float64,4)
fine_reffe   = ReferenceFE(QUAD,lagrangian,Float64,2)

V = FESpace(grid,fine_reffe)
cell_dof_ids = get_cell_dof_ids(V)
fine_face_nodes = fine_reffe.face_nodes

ftopo = get_grid_topology(grid)
fine_c2e_map = get_faces(ftopo,2,1)
fine_e2c_map = get_faces(ftopo,1,2)

# [Coarse edge][Fine child edge][Owned fine dofs]
c2f_edge_to_fine_dof = Vector{Vector{Vector{Int32}}}(undef,4)
for cE in 1:4 # For each coarse edge
  fedge_to_fine_dof = Vector{Vector{Int32}}(undef,length(c2f_edges[cE]))
  for (child_id,fE) in enumerate(c2f_edges[cE]) # For each child (fine edge)
    owner_cell = first(fine_e2c_map[fE])                      # Fine cell owner of child
    local_id   = findfirst(E->E==fE,fine_c2e_map[owner_cell]) # Local id (within owner cell) of the child
    fine_dofs  = cell_dof_ids[owner_cell][fine_face_nodes[4+local_id]]
    fedge_to_fine_dof[child_id] = fine_dofs
  end
  c2f_edge_to_fine_dof[cE] = fedge_to_fine_dof
end

coarse_face_nodes = coarse_reffe.face_nodes
c_edge_to_coarse_dof = coarse_face_nodes[5:8]

for cE in 1:4 # For each coarse edge
  coarse_dofs = c_edge_to_coarse_dof[cE]
  fine_dofs   = c2f_edge_to_fine_dof[cE]
  println(coarse_dofs, " <==> ",fine_dofs)
end

coarse_terms = _get_terms(SEGMENT,[4])
fine_terms = _get_terms(SEGMENT,[2])

coarse_dofs_above_fine_dofs = copy(c2f_edge_to_fine_dof)
for cE in 1:4
  coarse_dofs = c_edge_to_coarse_dof[cE]
  fine_dofs   = coarse_dofs_above_fine_dofs[cE]

  reordered_coarse_dofs = copy(coarse_dofs)
  reordered_coarse_dofs[coarse_terms] = coarse_dofs
  println(coarse_dofs)
  println(reordered_coarse_dofs)
  println(repeat('-',30))

  fine_dofs[1][fine_terms] = reordered_coarse_dofs[1:3]
  fine_dofs[2][fine_terms] = reordered_coarse_dofs[3:5]
end



end