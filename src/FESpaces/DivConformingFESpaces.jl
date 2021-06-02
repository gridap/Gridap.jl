# This source file is though to put that code required in order to
# customize ConformingFESpaces.jl to H(div)-conforming global FE Spaces built
# out of RaviartThomas FEs. In particular, this customization is in the
# definition of the shape functions (get_cell_shapefuns) and the DoFs
# (get_cell_dof_basis) of the **global** FE space, which requires a sign flip
# for those sitting on facets of the slave cell of the facet.

# Two key ingredients in the implementation of this type of ReferenceFE are the
# get_cell_shapefuns(model,cell_reffes,::Conformity) and
# get_cell_dof_basis(mode,cell_reffes,::Conformity) overloads.
# These are written such that, for each cell K, they return the shape functions
# and dof values in the *global* RT space. For a dof owned by a face which is shared by
# two cells, there is a master and a slave cell. The slave cell first computes the
# shape functions and dof values using local-to-cell data structures, but then flips the
# sign of both in order to get their corresponding counterparts in the **global**
# RT space. As as result we have the following:

# * When we interpolate a function into the global FE space, and we perform the cell-wise
#   DoF values to global DoF values gather operation, we can either extract the global DoF value
#   from the master or slave cell without worrying about the sign.
# * When we evaluate a global FE function, and we perform the global DoF values to
#   cell-wise DoF values scatter operation, we don't have to worry about the sign either.
#   On the slave cell, we will have both the sign of the DoF value, and the sign of the
#   shape function corresponding to the global DoF.
# * We do NOT have to use the signed determinant, but its absolute value, in the Piola Map.

function get_cell_dof_basis(model::DiscreteModel,
                       cell_reffe::AbstractArray{<:GenericRefFE{RaviartThomas}},
                       ::DivConformity)
    sign_flip = get_sign_flip(model, cell_reffe)
    lazy_map(_transform_rt_dof_basis,
             cell_reffe,
             get_cell_map(Triangulation(model)),
             sign_flip)
end

function get_cell_shapefuns(model::DiscreteModel,
                            cell_reffe::AbstractArray{<:GenericRefFE{RaviartThomas}},
                            ::DivConformity)
    sign_flip = get_sign_flip(model, cell_reffe)
    lazy_map(_transform_rt_shapefuns,
             cell_reffe,
             get_cell_map(Triangulation(model)),
             lazy_map(Broadcasting(constant_field), sign_flip))
end

struct SignFlipMap{T} <: Map
  model::T
end

function return_cache(k::SignFlipMap,reffe,cell_id)
  model = k.model
  D = num_cell_dims(model)
  gtopo = get_grid_topology(model)

  # Extract composition among cells and facets
  cell_wise_facets_ids = get_faces(gtopo, D, D - 1)
  cache_cell_wise_facets_ids = array_cache(cell_wise_facets_ids)

  # Extract cells around facets
  cells_around_facets = get_faces(gtopo, D - 1, D)
  cache_cells_around_facets = array_cache(cells_around_facets)

  (cell_wise_facets_ids,
   cache_cell_wise_facets_ids,
   cells_around_facets,
   cache_cells_around_facets,
   CachedVector(Bool))

end

function evaluate!(cache,k::SignFlipMap,reffe,cell_id)
  model = k.model

  cell_wise_facets_ids,
  cache_cell_wise_facets_ids,
  cells_around_facets,
  cache_cells_around_facets,
  sign_flip_cached = cache

  setsize!(sign_flip_cached, (num_dofs(reffe),))
  sign_flip = sign_flip_cached.array
  sign_flip .= false

  D = num_dims(reffe)
  face_own_dofs = get_face_own_dofs(reffe)
  facet_lid = get_offsets(get_polytope(reffe))[D] + 1
  cell_facets_ids = getindex!(cache_cell_wise_facets_ids,
                              cell_wise_facets_ids,
                              cell_id)
  for facet_gid in cell_facets_ids
      facet_cells_around = getindex!(cache_cells_around_facets,
           cells_around_facets,
           facet_gid)
      is_slave = (findfirst((x) -> (x == cell_id), facet_cells_around) == 2)
      if is_slave
          for dof in face_own_dofs[facet_lid]
              sign_flip[dof] = true
          end
      end
      facet_lid = facet_lid + 1
  end
  sign_flip
end

function get_sign_flip(model::DiscreteModel,
                       cell_reffe::AbstractArray{<:GenericRefFE{RaviartThomas}})
    lazy_map(SignFlipMap(model),
            cell_reffe,
            get_cell_to_bgcell(model))
end

function _transform_rt_shapefuns(reffe::GenericRefFE{RaviartThomas},
                                 phi::Field,
                                 sign_flip::AbstractVector{<:Field})
  ContraVariantPiolaMap()(get_shapefuns(reffe),phi,sign_flip)
end

function _transform_rt_dof_basis(reffe::GenericRefFE{RaviartThomas},
                       phi::Field,
                       sign_flip::AbstractVector{Bool})
  cache = return_cache(_transform_rt_dof_basis,reffe,phi,sign_flip)
  evaluate!(cache,_transform_rt_dof_basis,reffe,phi,sign_flip)
end

function return_cache(::typeof(_transform_rt_dof_basis),
                      reffe::GenericRefFE{RaviartThomas},
                      phi::Field,
                      ::AbstractVector{Bool})
  p = get_polytope(reffe)
  prebasis = get_prebasis(reffe)
  order = get_order(prebasis)
  et = return_type(prebasis)
  dofs = get_dof_basis(reffe)
  nodes, nf_nodes, nf_moments =  get_nodes(dofs),
                                 get_face_nodes_dofs(dofs),
                                 get_face_moments(dofs)
  db = MomentBasedDofBasis(nodes,nf_moments,nf_nodes)

  # Determine the type of the elements
  # in the range of Field. This is the
  # element type of the face_moments array.
  D = num_dims(reffe)
  pt = Point{D,et}
  T=return_type(phi,zero(pt))
  face_moments = [ similar(i,T)  for i in nf_moments ]

  Jt_q_cache = return_cache(∇(phi),db.nodes)
  cache = (db.nodes, db.face_nodes, nf_moments, face_moments, Jt_q_cache)
  cache
end


function evaluate!(cache,
                   ::typeof(_transform_rt_dof_basis),
                   reffe::GenericRefFE{RaviartThomas},
                   phi::Field,
                   sign_flip::AbstractVector{Bool})
  nodes, nf_nodes, nf_moments, face_moments, Jt_q_cache = cache

  Jt_q = evaluate!(Jt_q_cache,∇(phi),nodes)
  face_own_dofs=get_face_own_dofs(reffe)
  for face in 1:length(face_moments)
    moments = nf_moments[face]
    if length(moments) > 0
      sign = (-1)^sign_flip[face_own_dofs[face][1]]
      num_qpoints, num_moments = size(moments)
      for i in 1:num_qpoints
        Jt_q_i = Jt_q[nf_nodes[face][i]]
        change = sign * meas(Jt_q_i) * pinvJt(Jt_q_i)
        for j in 1:num_moments
          face_moments[face][i,j] = change ⋅ nf_moments[face][i,j]
        end
      end
    end
  end
  MomentBasedDofBasis(nodes,face_moments,nf_nodes)
end
