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

struct TransformRTDofBasis{Dc,Dp} <: Map end

function get_cell_dof_basis(
  model::DiscreteModel{Dc,Dp},
  cell_reffe::AbstractArray{<:GenericRefFE{<:DivConforming}},
  ::DivConformity,
  sign_flip = get_sign_flip(model, cell_reffe)
) where {Dc,Dp}
  cell_map  = get_cell_map(get_grid(model))
  Jt  = lazy_map(Broadcasting(∇),cell_map)
  x   = lazy_map(get_nodes,lazy_map(get_dof_basis,cell_reffe))
  Jtx = lazy_map(evaluate,Jt,x)
  k = TransformRTDofBasis{Dc,Dp}()
  lazy_map(k,cell_reffe,Jtx,sign_flip)
end

function get_cell_shapefuns(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{<:DivConforming}},
  ::DivConformity,
  sign_flip = get_sign_flip(model, cell_reffe)
)
  cell_map = get_cell_map(get_grid(model))
  cell_shapefuns = lazy_map(get_shapefuns,cell_reffe)
  k = ContraVariantPiolaMap()
  lazy_map(k,cell_shapefuns,cell_map,lazy_map(Broadcasting(constant_field),sign_flip))
end

struct SignFlipMap{T} <: Map
  model::T
  facet_owners::Vector{Int32}
end

function SignFlipMap(model)
  facet_owners = compute_facet_owners(model)
  SignFlipMap(model,facet_owners)
end

function return_cache(k::SignFlipMap,reffe,facet_own_dofs,cell)
  model = k.model
  Dc = num_cell_dims(model)
  topo = get_grid_topology(model)

  cell_facets = get_faces(topo, Dc, Dc-1)
  cell_facets_cache = array_cache(cell_facets)

  return cell_facets,cell_facets_cache,CachedVector(Bool)
end

function evaluate!(cache,k::SignFlipMap,reffe,facet_own_dofs,cell)
  cell_facets,cell_facets_cache,sign_flip_cache = cache
  facet_owners = k.facet_owners

  setsize!(sign_flip_cache, (num_dofs(reffe),))
  sign_flip = sign_flip_cache.array
  sign_flip .= false

  facets = getindex!(cell_facets_cache,cell_facets,cell)
  for (lfacet,facet) in enumerate(facets)
    owner = facet_owners[facet]
    if owner != cell
      for dof in facet_own_dofs[lfacet]
        sign_flip[dof] = true
      end
    end
  end

  return sign_flip
end

function get_sign_flip(
  model::DiscreteModel{Dc},
  cell_reffe::AbstractArray{<:GenericRefFE{<:DivConforming}}
) where Dc
  # Comment: lazy_maps on cell_reffes are very optimised, since they are CompressedArray/FillArray
  get_facet_own_dofs(reffe) = view(get_face_own_dofs(reffe),get_dimrange(get_polytope(reffe),Dc-1))
  cell_facet_own_dofs = lazy_map(get_facet_own_dofs,cell_reffe)
  cell_ids = IdentityVector(Int32(num_cells(model)))
  lazy_map(SignFlipMap(model),cell_reffe,cell_facet_own_dofs,cell_ids)
end

function compute_facet_owners(model::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  topo = get_grid_topology(model)
  facet_to_cell = get_faces(topo, Dc-1, Dc)

  nfacets = num_faces(topo, Dc-1)
  owners = Vector{Int32}(undef, nfacets)
  for facet in 1:nfacets
    facet_cells = view(facet_to_cell, facet)
    owners[facet] = first(facet_cells)
  end

  return owners
end

function return_cache(
  ::TransformRTDofBasis{Dc,Dp},
  reffe::GenericRefFE{<:DivConforming},
  Jtx,
  ::AbstractVector{Bool}
) where {Dc,Dp}
  # @santiagobadia: Hack as above
  et = eltype(return_type(get_prebasis(reffe)))
  dofs = get_dof_basis(reffe)

  nodes = get_nodes(dofs)
  nf_nodes = get_face_nodes_dofs(dofs)
  nf_moments = get_face_moments(dofs)
  db = MomentBasedDofBasis(nodes,nf_moments,nf_nodes)
  face_moments = [ similar(i,VectorValue{Dp,et})  for i in nf_moments ]

  cache = (db.nodes, db.face_nodes, nf_moments, face_moments)
  cache
end

function evaluate!(
  cache,
  ::TransformRTDofBasis,
  reffe::GenericRefFE{<:DivConforming},
  Jt_q,
  sign_flip::AbstractVector{Bool}
)
  nodes, nf_nodes, nf_moments, face_moments = cache
  face_own_dofs=get_face_own_dofs(reffe)
  for face in 1:length(face_moments)
    nf_moments_face   = nf_moments[face]
    face_moments_face = face_moments[face]
    if length(nf_moments_face) > 0
      sign = (-1)^sign_flip[face_own_dofs[face][1]]
      num_qpoints, num_moments = size(nf_moments_face)
      for i in 1:num_qpoints
        Jt_q_i = Jt_q[nf_nodes[face][i]]
        change = sign * meas(Jt_q_i) * pinvJt(Jt_q_i)
        for j in 1:num_moments
          face_moments_face[i,j] = change ⋅ nf_moments_face[i,j]
        end
      end
    end
  end
  MomentBasedDofBasis(nodes,face_moments,nf_nodes)
end


# Support for DIV operator
function DIV(f::LazyArray{<:Fill{T}}) where T
  df=DIV(f.args[1])
  k=f.maps.value
  lazy_map(k,df)
end
function DIV(f::LazyArray{<:Fill{Broadcasting{Operation{ContraVariantPiolaMap}}}})
  ϕrgₖ       = f.args[1]
  fsign_flip = f.args[4]
  div_ϕrgₖ = lazy_map(Broadcasting(divergence),ϕrgₖ)
  fsign_flip=lazy_map(Broadcasting(Operation(x->(-1)^x)), fsign_flip)
  lazy_map(Broadcasting(Operation(*)),fsign_flip,div_ϕrgₖ)
end
function DIV(a::LazyArray{<:Fill{typeof(linear_combination)}})
  i_to_basis = DIV(a.args[2])
  i_to_values = a.args[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end
