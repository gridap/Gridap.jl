
function get_cell_dof_basis(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}},
  ::CurlConformity)

  Dc        = num_dims(reffe)
  Dp        = first(size(return_type(phi,zero(pt))))
  cell_dofs = lazy_map(get_dof_basis,cell_reffe)
  cell_ownids = lazy_map(get_face_own_dofs,cell_reffe)
  cell_map = get_cell_map(Triangulation(model))
  cell_Jt = lazy_map(Broadcasting(∇),cell_map)
  cell_x = lazy_map(get_nodes,cell_dofs)
  cell_Jtx  = lazy_map(evaluate,cell_Jt,cell_x)
  k = TransformNedelecDofBasis{Dc,Dp}()
  lazy_map(k,cell_dofs,cell_Jtx,cell_ownids)
end

struct TransformNedelecDofBasis{Dc,Dp} <: Map end;

  function return_cache(::TransformNedelecDofBasis{Dc,Dp},dofs,Jtx,ownids) where {Dc,Dp}
    # Assumes the same element type in all the mesh
    nodes, nf_nodes, nf_moments = get_nodes(dofs),get_face_nodes_dofs(dofs),get_face_moments(dofs)
    face_moments = [ similar(i,VectorValue{Dp,typeof(dofs.nodes[1][1])})  for i in nf_moments ]
    db = MomentBasedDofBasis(nodes,face_moments,nf_nodes)
    db
  end

function evaluate!(cache,::TransformNedelecDofBasis,dofs,Jtx,ownids)
  basis = cache # Assumes the same element type in all the mesh
  for face in 1:length(ownids)
    face_dofs_ids = ownids[face]
    face_point_ids = dofs.face_nodes[face]
    face_moments_out = basis.face_moments[face]
    face_moments_in = dofs.face_moments[face]
    for p in 1:length(face_point_ids)
      F = transpose(Jtx[p])
      for i in 1:length(face_dofs_ids)
        face_moments_out[p,i] = F⋅face_moments_in[p,i]
      end
    end
  end
  basis
end

function get_cell_shapefuns(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}},
  ::CurlConformity)

  cell_reffe_shapefuns = lazy_map(get_shapefuns,cell_reffe)
  cell_map = get_cell_map(Triangulation(model))
  cell_Jt  = lazy_map(Broadcasting(∇),cell_map)
  k = ReferenceFEs.CoVariantPiolaMap()
  lazy_map(k,cell_reffe_shapefuns,cell_Jt)
end
