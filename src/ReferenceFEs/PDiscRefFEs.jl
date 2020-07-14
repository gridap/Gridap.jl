
function _PDiscRefFE(::Type{T},p::Polytope,orders) where T
    order = first(orders)
    @notimplementedif any( orders .!= order ) "Anisotropic serentopity FEs not allowed"
    _PDiscRefFE(T,p,order)
end

function _PDiscRefFE(::Type{T},p::Polytope,order::Integer) where T
  D = num_cell_dims(p)
  extrusion = tfill(TET_AXIS,Val{D}())
  simplex = ExtrusionPolytope(extrusion)
  reffe = LagrangianRefFE(T,simplex,order)
  metadata = nothing
  face_nodes = [Int[] for face in 1:num_faces(p)]
  face_nodes[end] = collect(1:num_nodes(reffe))
  dofs = get_dof_basis(reffe)
  face_dofs = _generate_face_own_dofs(face_nodes,dofs.node_and_comp_to_dof)

  reffe = GenericRefFE(
    num_dofs(reffe),
    p,
    get_prebasis(reffe),
    get_dof_basis(reffe),
    L2Conformity(),
    metadata,
    face_dofs)
  GenericLagrangianRefFE(reffe,face_nodes)
end

