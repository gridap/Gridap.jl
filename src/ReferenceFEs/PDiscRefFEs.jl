"""
    struct PDiscRefFE{D} <: NodalReferenceFE{D}
      # Private fields
    end
"""
struct PDiscRefFE{D} <: NodalReferenceFE{D}
  p_reffe::LagrangianRefFE{D}
  polytope::Polytope{D}
  function PDiscRefFE(::Type{T},p::Polytope,order::Integer) where T
    D = num_cell_dims(p)
    extrusion = tfill(TET_AXIS,Val{D}())
    simplex = ExtrusionPolytope(extrusion)
    p_reffe = LagrangianRefFE(T,simplex,order)
    new{D}(p_reffe,p)
  end
end

# ReferenceFE interface

num_dofs(reffe::PDiscRefFE) = num_dofs(reffe.p_reffe)

get_polytope(reffe::PDiscRefFE) = reffe.polytope

get_prebasis(reffe::PDiscRefFE) = get_prebasis(reffe.p_reffe)

get_dof_basis(reffe::PDiscRefFE) = get_dof_basis(reffe.p_reffe)

function get_face_own_dofs(reffe::PDiscRefFE)
  _compute_non_conf_face_own_x(reffe,num_dofs(reffe))
end

function get_face_own_dofs_permutations(reffe::PDiscRefFE)
  _compute_non_conf_face_own_x_permutations(reffe,num_dofs(reffe))
end

get_face_dofs(reffe::PDiscRefFE) = get_face_own_dofs(reffe)

# Nodal ReferenceFE interface

get_node_coordinates(reffe::PDiscRefFE) = get_node_coordinates(reffe.p_reffe)

get_node_and_comp_to_dof(reffe::PDiscRefFE) = get_node_and_comp_to_dof(reffe.p_reffe)

function get_face_own_nodes(reffe::PDiscRefFE)
  _compute_non_conf_face_own_x(reffe,num_nodes(reffe))
end

function get_face_own_nodes_permutations(reffe::PDiscRefFE)
  _compute_non_conf_face_own_x_permutations(reffe,num_nodes(reffe))
end

get_face_nodes(reffe::PDiscRefFE) = get_face_own_nodes(reffe)

# Some more API for this reffe type

get_dof_to_node(reffe::PDiscRefFE) = get_dof_to_node(reffe.p_reffe)

get_dof_to_comp(reffe::PDiscRefFE) = get_dof_to_comp(reffe.p_reffe)

# Helpers

function _compute_non_conf_face_own_x(reffe,ndofs)
  polytope = get_polytope(reffe)
  nfaces = num_faces(polytope)
  dofs = [ Int[]  for face in 1:nfaces]
  dofs[end] = collect(1:ndofs)
  dofs
end

function _compute_non_conf_face_own_x_permutations(reffe,ndofs)
  polytope = get_polytope(reffe)
  nfaces = num_faces(polytope)
  face_to_perms = get_face_vertex_permutations(polytope)
  face_to_pindex_to_dofs = Vector{Vector{Int}}[]
  for face in 1:nfaces
    nperms = length(face_to_perms[face])
    if face != nfaces
      pindex_to_dofs = [Int[] for i in 1:nperms]
    else
      pindex_to_dofs = [fill(INVALID_PERM,ndofs) for i in 1:nperms]
      pindex_to_dofs[1] = collect(1:ndofs)
    end
    push!(face_to_pindex_to_dofs,pindex_to_dofs)
  end
  face_to_pindex_to_dofs
end

