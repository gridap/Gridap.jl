
function _PDiscRefFE(::Type{V},p::Polytope,orders, poly_type) where V
    order = first(orders)
    @notimplementedif any( orders .!= order ) "Anisotropic serendipity FEs not allowed"
    _PDiscRefFE(V,p,order, poly_type)
end

function _PDiscRefFE(::Type{V},p::Polytope,order::Integer, poly_type) where V

  D = num_cell_dims(p)
  extrusion = tfill(TET_AXIS,Val{D}())
  simplex = ExtrusionPolytope(extrusion)
  reffe = LagrangianRefFE(V,simplex,order)
  metadata = nothing
  face_nodes = [Int[] for face in 1:num_faces(p)]
  face_nodes[end] = collect(1:num_nodes(reffe))
  # basis for ( SᵣΛᴰ(□ᴰ) )ⁿ, n=num_indep_components(V)
  basis = FEEC_poly_basis(Val(D),V,order,D,:S,poly_type; cart_prod=(V <: MultiValue))
  dofs = get_dof_basis(reffe)
  face_own_dofs = _generate_face_own_dofs(face_nodes,dofs.node_and_comp_to_dof)

  conf = L2Conformity()

  reffe = GenericRefFE{typeof(conf)}(
    num_dofs(reffe),
    p,
    basis, # pre-basis
    dofs,
    conf,
    metadata,
    face_own_dofs)
  GenericLagrangianRefFE(reffe,face_nodes)
end

