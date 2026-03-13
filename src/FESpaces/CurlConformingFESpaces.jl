

function _no_sign_flip(model, cell_reffe)
    # Assuming that the first element exists, 
    # and that all the elements are of the same type
    cell_reffe_one=first(cell_reffe)
    sign_flip_single_cell = Falses(num_dofs(cell_reffe_one))
    Fill(sign_flip_single_cell, num_cells(model))
end


# No sign flip to be applied in 2D for Nedelec facet DoFs
function get_sign_flip(model::DiscreteModel{2,Dp},
                       cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}}) where Dp
    _no_sign_flip(model, cell_reffe)
end


function get_sign_flip(model::DiscreteModel{3,Dp},
                       cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}}) where Dp                   
    p=get_polytope(first(cell_reffe))
    if (is_simplex(p))
        # No sign flip to be applied to tetrahedral meshes for Nedelec facet DoFs
        _no_sign_flip(model, cell_reffe)
    else
       @assert is_n_cube(p)
       lazy_map(SignFlipMap(model),
                cell_reffe,
                IdentityVector(Int32(num_cells(model))))
    end
end


function get_cell_dof_basis(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}},
  ::CurlConformity,
  sign_flip=get_sign_flip(model, cell_reffe))
  cell_map  = get_cell_map(Triangulation(model))
  phi       = cell_map[1]
  reffe     = cell_reffe[1]
  Dc        = num_dims(reffe)
  et        = eltype(return_type(get_prebasis(reffe)))
  pt        = Point{Dc,et}
  Dp        = first(size(return_type(phi,zero(pt))))
  cell_dofs = lazy_map(get_dof_basis,cell_reffe)
  cell_ownids = lazy_map(get_face_own_dofs,cell_reffe)
  cell_map = get_cell_map(Triangulation(model))
  cell_Jt = lazy_map(Broadcasting(∇),cell_map)
  cell_x = lazy_map(get_nodes,cell_dofs)
  cell_Jtx  = lazy_map(evaluate,cell_Jt,cell_x)
  k = TransformNedelecDofBasis{Dc,Dp}()
  lazy_map(k,cell_dofs,cell_Jtx,cell_ownids,sign_flip)
end

struct TransformNedelecDofBasis{Dc,Dp} <: Map end;

  function return_cache(::TransformNedelecDofBasis{Dc,Dp},
                        dofs,
                        Jtx,
                        ownids,
                        ::AbstractVector{Bool}) where {Dc,Dp}
    # Assumes the same element type in all the mesh
    nodes, nf_nodes, nf_moments = get_nodes(dofs),get_face_nodes_dofs(dofs),get_face_moments(dofs)
    face_moments = [ similar(i,VectorValue{Dp,typeof(dofs.nodes[1][1])})  for i in nf_moments ]
    db = MomentBasedDofBasis(nodes,face_moments,nf_nodes)
    db
  end

function evaluate!(cache,
                   ::TransformNedelecDofBasis{Dc,Dp},
                   dofs,
                   Jtx,
                   ownids,
                   sign_flip::AbstractVector{Bool}) where {Dc,Dp}
  basis = cache # Assumes the same element type in all the mesh

  for face in 1:length(ownids)
    face_dofs_ids = ownids[face]
    face_point_ids = dofs.face_nodes[face]
    face_moments_out = basis.face_moments[face]
    face_moments_in = dofs.face_moments[face]
    for p in 1:length(face_point_ids)
      F = transpose(Jtx[face_point_ids[p]])
      for i in 1:length(face_dofs_ids)
        sign = (-1)^sign_flip[face_dofs_ids[i]]
        face_moments_out[p,i] = sign * (F⋅face_moments_in[p,i])
      end
    end
  end
  basis
end

function get_cell_shapefuns(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{Nedelec}},
  ::CurlConformity,
  sign_flip=get_sign_flip(model, cell_reffe))

  cell_reffe_shapefuns = lazy_map(get_shapefuns,cell_reffe)
  cell_map = get_cell_map(Triangulation(model))
  k = CoVariantPiolaMap()
  lazy_map(k,
           cell_reffe_shapefuns,
           cell_map,
           lazy_map(Broadcasting(constant_field), sign_flip))
end

