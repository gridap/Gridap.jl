
struct FineToCoarseIndices <: AbstractVector{Tuple{Vector{Int32},Vector{Int32}}}
  fcell_to_cids :: Vector{Vector{Int32}}
  cid_to_fcells :: Vector{Vector{Int32}}
  cid_to_fids   :: Vector{Vector{Int32}}
  one_to_one    :: Bool
  function FineToCoarseIndices(
    fcell_to_cids::AbstractVector{<:AbstractVector{<:Integer}}
  )
    fcell_to_cids = collect(Vector{Int32},fcell_to_cids)

    n_cids = maximum(map(maximum,fcell_to_cids))
    cid_to_fcells = [Int32[] for i in 1:n_cids]
    cid_to_fids = [Int32[] for i in 1:n_cids]
    for (fcell,cids) in enumerate(fcell_to_cids)
      for (fid,cid) in enumerate(cids)
        push!(cid_to_fcells[cid],fcell)
        push!(cid_to_fids[cid],fid)
      end
    end
    one_to_one = all(fids -> isone(length(fids)),cid_to_fids)
    new(fcell_to_cids,cid_to_fcells,cid_to_fids,one_to_one)
  end
end

Base.size(a::FineToCoarseIndices) = (length(a.cid_to_fids),)
function Base.getindex(a::FineToCoarseIndices,cid::Integer)
  (a.cid_to_fcells[cid],a.cid_to_fids[cid])
end

# Note for developpers:
# We use FineToCoarseArray to represent many data-structures, from quadrature points to 
# finite-element basis. For performance reasons, we sometimes want to keep a copy of the 
# coarse data (but not always!). For instance: 
#  - for quad points, we want to keep the coarse points. Otherwise, we woudl have to 
#    recompute them for every cell when the CompositeQuadrature is used with other 
#    gridap basis. 
#  - for finite-element basis, we do NOT want to generate the coarse basis, otherwise 
#    we would be creating it many times (when it should never be used). In this case, we 
#    create it on the go (see getindex! specialisations).
# All in all, this is why coarse_data can be of type Nothing.
struct FineToCoarseArray{T,A,B,C} <: AbstractVector{T}
  rrule       :: A
  coarse_data :: B
  fine_data   :: C
  ids         :: FineToCoarseIndices

  function FineToCoarseArray{T}(
    rrule::RefinementRule,
    coarse_data::Union{Vector{T},Nothing},
    fine_data::AbstractVector{<:AbstractVector},
    ids::FineToCoarseIndices
  ) where T
    A = typeof(rrule)
    B = typeof(coarse_data)
    C = typeof(fine_data)
    new{T,A,B,C}(rrule,coarse_data,fine_data,ids)
  end
end

function FineToCoarseArray(
  rrule::RefinementRule,
  fine_data::AbstractVector{<:AbstractVector},
  ids::FineToCoarseIndices
)
  T = combine_fine_to_coarse_type(rrule,fine_data,ids)
  return FineToCoarseArray{T}(rrule,nothing,fine_data,ids)
end

function FineToCoarseArray(
  rrule::RefinementRule,
  fine_data::AbstractVector{<:AbstractVector},
  connectivity::AbstractVector{<:AbstractVector{<:Integer}}
)
  ids = FineToCoarseIndices(connectivity)
  return FineToCoarseArray(rrule,fine_data,ids)
end

function FineToCoarseArray(
  rrule::RefinementRule,
  fine_data::AbstractVector{<:AbstractVector}
)
  offsets = cumsum(map(length,fine_data)) .- length(first(fine_data)) .+ 1
  connectivity = map(fine_data,offsets) do fdata,o
    collect(Int32,o:o+length(fdata)-1)
  end
  return FineToCoarseArray(rrule,fine_data,connectivity)
end

Base.size(a::FineToCoarseArray) = size(a.ids)

Base.getindex(a::FineToCoarseArray{T},i::Integer) where T = getindex(a.coarse_data,i)

function Base.getindex(a::FineToCoarseArray{T,A,Nothing},i::Integer) where {T,A}
  fcells, fids = getindex(a.ids,i)
  fdata = map((fcell,fid) -> getindex(a.fine_data[fcell],Int(fid)),fcells,fids)
  return combine_fine_to_coarse(a.rrule,fdata,fcells)
end

function Arrays.get_children(n::TreeNode,a::FineToCoarseArray)
  (similar_tree_node(n,a.rrule),similar_tree_node(n,a.fine_data))
end

@inline function combine_fine_to_coarse_type(
  rr::RefinementRule,fine_data::AbstractVector{<:AbstractVector{T}},ids::FineToCoarseIndices
) where T <: Field
  return FineToCoarseField{Vector{T},typeof(rr),Vector{Int32}}
end
@inline function combine_fine_to_coarse(
  rr::RefinementRule,fine_fields::Vector{<:Field},child_ids::Vector{<:Integer}
)
  FineToCoarseField(fine_fields,rr,child_ids)
end

struct FineToCoarseDof <: Dof end # Should we implement this properly?
@inline function combine_fine_to_coarse_type(
  rr::RefinementRule,fine_data::AbstractVector{<:AbstractVector{T}},ids::FineToCoarseIndices
) where T <: Dof
  return FineToCoarseDof
end
@inline function combine_fine_to_coarse(
  rr::RefinementRule,fine_dofs::Vector{<:Dof},child_ids::Vector{<:Integer}
)
  FineToCoarseDof()
end

@inline function combine_fine_to_coarse_type(
  rr::RefinementRule,fine_data::AbstractVector{<:AbstractVector{T}},ids::FineToCoarseIndices
) where T <: Point
  return T
end
@inline function combine_fine_to_coarse(
  rr::RefinementRule,fine_pts::Vector{<:Point},child_ids::Vector{<:Integer}
)
  cmaps = get_cell_map(rr)
  evaluate(cmaps[first(child_ids)],first(fine_pts))
end

@inline function combine_fine_to_coarse_type(
  rr::RefinementRule,fine_data::AbstractVector{<:AbstractVector{T}},ids::FineToCoarseIndices
) where T <: Real
  return T
end
@inline function combine_fine_to_coarse(
  rr::RefinementRule,fine_vals::Vector{<:Real},child_ids::Vector{<:Integer}
)
  sum(fine_vals)
end

# MacroDofBasis evaluation

const MacroDofBasis = FineToCoarseArray{<:Dof}
const MacroFEBasis = FineToCoarseArray{<:Field}

function Arrays.return_type(a::MacroDofBasis,b::MacroFEBasis)
  @check a.rrule == b.rrule
  types = map(return_type,a.fine_data,b.fine_data)
  T = typejoin(types...)
  @assert all(t -> t <: T, types)
  return T
end

function Arrays.return_value(a::MacroDofBasis,b::MacroFEBasis)
  @check a.rrule == b.rrule
  T = return_type(a,b)
  return similar(T,length(a),length(b))
end

function Arrays.return_cache(a::MacroDofBasis,b::MacroFEBasis)
  @check a.rrule == b.rrule
  caches = map(return_cache,a.fine_data,b.fine_data)
  
  T = eltype(evaluate!(first(caches),first(a.fine_data),first(b.fine_data)))
  res = zeros(T,length(a),length(b))
  return res, caches
end

function Arrays.evaluate!(cache,a::MacroDofBasis,b::MacroFEBasis)
  res, caches = cache
  fill!(res,zero(eltype(res)))
  cell_vals = map(evaluate!,caches,a.fine_data,b.fine_data)
  for fcell in 1:num_subcells(a.rrule)
    I = a.ids.fcell_to_cids[fcell]
    J = b.ids.fcell_to_cids[fcell]
    res[I,J] .= cell_vals[fcell]
  end
  return res
end

function Arrays.return_value(dofs::MacroDofBasis,f::Field)
  cmap = get_cell_map(dofs.rrule)
  v = return_value(first(dofs.fine_data),f∘first(cmap))
  return similar(v,length(dofs))
end

function Arrays.return_cache(dofs::MacroDofBasis,f::Field)
  cmap = get_cell_map(dofs.rrule)
  caches = map((dofs_k,mk) -> return_cache(dofs_k,f∘mk),dofs.fine_data,cmap)
  T = eltype(evaluate!(first(caches),first(dofs.fine_data),f∘first(cmap)))
  res = zeros(T,length(dofs))
  return res, cmap, caches
end

function Arrays.evaluate!(cache,dofs::MacroDofBasis,f::Field)
  res, cmap, caches = cache
  fill!(res,zero(eltype(res)))
  cell_vals = map((cache_k,dofs_k,mk) -> evaluate!(cache_k,dofs_k,f∘mk),caches,dofs.fine_data,cmap)
  for fcell in 1:num_subcells(dofs.rrule)
    dof_ids = dofs.ids.fcell_to_cids[fcell]
    res[dof_ids] .= cell_vals[fcell]
  end
  return res
end

# MacroFEBasis evaluation with coarse points

function Arrays.return_type(a::MacroFEBasis,xc::AbstractArray{<:Point})
  return_type(first(a.fine_data),xc)
end

function Arrays.return_value(a::MacroFEBasis,xc::AbstractArray{<:Point})
  TM = return_type(a,xc)
  return similar(TM,length(xc),length(a))
end

function Arrays.return_cache(a::MacroFEBasis,xc::AbstractArray{<:Point})
  rr = a.rrule

  k = CoarseToFinePointMap()
  geo_cache = return_cache(k,rr,xc)
  xf, ids = evaluate!(geo_cache,k,rr,xc)

  eval_caches = map(return_cache,a.fine_data,xf)

  T = eltype(evaluate!(first(eval_caches),first(a.fine_data),first(xf)))
  res_cache = CachedArray(zeros(T,length(xc),length(a)))
  return res_cache, k, geo_cache, eval_caches
end

function Arrays.evaluate!(caches, a::MacroFEBasis,xc::AbstractArray{<:Point})
  res_cache, k, geo_cache, eval_caches = caches
  setsize!(res_cache,(length(xc),length(a)))
  res = res_cache.array

  fill!(res,zero(eltype(res)))
  xf, ids = evaluate!(geo_cache,k,a.rrule,xc)
  
  for fcell in 1:num_subcells(a.rrule)
    r = xf.ptrs[fcell]:xf.ptrs[fcell+1]-1
    vals = evaluate!(eval_caches[fcell],a.fine_data[fcell],view(xf.data,r))
    I  = view(ids.data,r)
    J = a.ids.fcell_to_cids[fcell]
    res[I,J] .= vals
  end
  return res
end

# MacroFEBasis evaluation with fine points

function Arrays.return_type(a::MacroFEBasis,b::FineToCoarseArray{<:Point})
  @check a.rrule == b.rrule
  return return_type(first(a.fine_data),first(b.fine_data))
end

function Arrays.return_value(a::MacroFEBasis,b::FineToCoarseArray{<:Point})
  MT = return_type(a,b)
  return similar(MT,length(b),length(a))
end

function Arrays.return_cache(a::MacroFEBasis,b::FineToCoarseArray{<:Point})
  @check a.rrule == b.rrule
  caches = map(return_cache,a.fine_data,b.fine_data)
  
  T = eltype(evaluate!(first(caches),first(a.fine_data),first(b.fine_data)))
  res_cache = CachedArray(zeros(T,length(b),length(a)))
  return res_cache, caches
end

function Arrays.evaluate!(cache,a::MacroFEBasis,b::FineToCoarseArray{<:Point})
  res_cache, caches = cache
  setsize!(res_cache,(length(b),length(a)))
  res = res_cache.array

  fill!(res,zero(eltype(res)))
  for fcell in 1:num_subcells(a.rrule)
    vals = evaluate!(caches[fcell],a.fine_data[fcell],b.fine_data[fcell])
    I = b.ids.fcell_to_cids[fcell]
    J = a.ids.fcell_to_cids[fcell]
    res[I,J] .= vals
  end
  return res
end

# MacroFEBasis optimisations

function Fields.linear_combination(a::AbstractVector{<:Number},b::MacroFEBasis)
  rrule, ids = b.rrule, b.ids
  fcoeffs = map(cids -> view(a,cids), ids.fcell_to_cids)
  ffields = map(linear_combination,fcoeffs,b.fine_data)
  return FineToCoarseField(ffields,rrule)
end

function Arrays.evaluate!(cache,k::Broadcasting{typeof(Fields.∇)},a::MacroFEBasis)
  rrule  = a.rrule
  cell_maps = get_cell_map(rrule)
  cell_grads = lazy_map(k,a.fine_data)
  fields = lazy_map(push_∇,cell_grads,cell_maps)
  return FineToCoarseArray(a.rrule,fields,a.ids)
end

function Arrays.evaluate!(cache,k::Broadcasting{typeof(Fields.∇∇)},a::MacroFEBasis)
  rrule  = a.rrule
  cell_maps = get_cell_map(rrule)
  cell_grads = lazy_map(k,a.fine_data)
  fields = lazy_map(push_∇∇,cell_grads,cell_maps)
  return FineToCoarseArray(a.rrule,fields,a.ids)
end

############################################################################################
# MacroReferenceFE

struct MacroRefFE <: ReferenceFEName end

"""
    MacroReferenceFE(rrule::RefinementRule,reffes::AbstractVector{<:ReferenceFE})

Constructs a ReferenceFE for a macro-element, given a RefinementRule and a set of
ReferenceFEs for the subcells.

For performance, these should be paired with CompositeQuadratures.
"""
function MacroReferenceFE(
  rrule::RefinementRule,
  reffes::AbstractVector{<:ReferenceFE};
  conformity = Conformity(first(reffes))
)
  @check length(reffes) == num_subcells(rrule)

  grid = rrule.ref_grid
  space = FESpace(grid,reffes;conformity=conformity)

  conn = get_cell_dof_ids(space)
  basis = FineToCoarseArray(rrule,map(get_shapefuns,reffes),conn)
  dofs  = FineToCoarseArray(rrule,map(get_dof_basis,reffes),conn)
  face_dofs = get_cface_to_dofs(rrule,space,reffes)
  face_own_dofs = get_cface_to_own_dofs(rrule,space,reffes)
  face_own_perms = get_cface_to_own_dof_permutations(rrule,space,reffes)

  ndofs = num_free_dofs(space)
  poly = get_polytope(rrule)
  # This is a hack to be able to compute the orders
  prebasis = FineToCoarseArray(rrule,collect(map(get_prebasis,reffes)))
  metadata = (rrule,conn,face_own_dofs,face_own_perms)

  return GenericRefFE{MacroRefFE}(
    ndofs,poly,prebasis,dofs,conformity,metadata,face_dofs,basis
  )
end

ReferenceFEs.get_order(reffe::GenericRefFE{MacroRefFE}) = maximum(get_orders(reffe))

function ReferenceFEs.get_orders(reffe::GenericRefFE{MacroRefFE})
  prebasis = get_prebasis(reffe)
  subcell_prebasis = prebasis.fine_data
  subcell_orders = map(get_orders,subcell_prebasis)
  return map(maximum,subcell_orders)
end

function ReferenceFEs.get_face_own_dofs(reffe::GenericRefFE{MacroRefFE}, conf::Conformity)
  @check conf == Conformity(reffe)
  rrule,conn,face_own_dofs,face_own_perms = ReferenceFEs.get_metadata(reffe)
  return face_own_dofs
end

function ReferenceFEs.get_face_own_dofs_permutations(reffe::GenericRefFE{MacroRefFE}, conf::Conformity)
  @check conf == Conformity(reffe)
  rrule,conn,face_own_dofs,face_own_perms = ReferenceFEs.get_metadata(reffe)
  return face_own_perms
end

function ReferenceFEs.get_face_own_dofs(reffe::GenericRefFE{MacroRefFE}, ::L2Conformity)
  return ReferenceFEs._get_face_own_dofs_l2(reffe)
end

function ReferenceFEs.get_face_own_dofs_permutations(reffe::GenericRefFE{MacroRefFE}, ::L2Conformity)
  face_own_dofs = ReferenceFEs.get_face_own_dofs(reffe,L2Conformity())
  return ReferenceFEs._trivial_face_own_dofs_permutations(face_own_dofs)
end

function ReferenceFEs.Conformity(reffe::GenericRefFE{MacroRefFE},sym::Symbol)
  if sym == :L2
    L2Conformity()
  else
    Conformity(reffe)
  end
end

"""
    get_cface_to_own_dofs(
      rrule::RefinementRule, 
      space::FESpace, 
      reffes::AbstractVector{<:ReferenceFE}
    )

Given a RefinementRule and a FESpace defined on it's grid (with a set of ReferenceFEs),
returns the dofs owned by each face of the macro-element (i.e each face of the
underlying polytope).
"""
function get_cface_to_own_dofs(
  rrule::RefinementRule{<:Polytope{Dc}},
  space::FESpace,
  reffes::AbstractVector{<:ReferenceFE}
) where Dc
  cface_to_own_fface_to_own_dofs = get_cface_to_own_fface_to_own_dofs(rrule,space,reffes)
  cface_to_own_dofs = map(x -> vcat(x...), cface_to_own_fface_to_own_dofs)
  return cface_to_own_dofs
end

"""
    get_cface_to_dofs(
      rrule::RefinementRule, 
      space::FESpace, 
      reffes::AbstractVector{<:ReferenceFE}
    )

Given a RefinementRule and a FESpace defined on it's grid (with a set of ReferenceFEs),
returns the dofs in the closure of each face of the macro-element (i.e each face of the
underlying polytope).
"""
function get_cface_to_dofs(
  rrule::RefinementRule{<:Polytope{Dc}},
  space::FESpace,
  reffes::AbstractVector{<:ReferenceFE}
) where Dc
  cface_to_own_fface_to_own_dofs = get_cface_to_own_fface_to_own_dofs(rrule,space,reffes)
  cface_to_fface_to_dofs = aggregate_cface_to_own_fface_data(
    rrule, cface_to_own_fface_to_own_dofs
  )
  cface_to_dofs = map(x -> vcat(x...), cface_to_fface_to_dofs)
  return cface_to_dofs
end

"""
    get_cface_to_own_fface_to_own_dofs(
      rrule::RefinementRule, space::FESpace, reffes::AbstractVector{<:ReferenceFE}
    )

Given a RefinementRule and a FESpace defined on it's grid (with a set of ReferenceFEs),
returns for each coarse face the dofs owned by each fine face owned by the coarse face.

Returns two arrays

- [coarse face][owned fine face] -> [owned dofs]
- [coarse face][owned fine face][local permutation] -> [local permuted dofs]

"""
function get_cface_to_own_fface_to_own_dofs(
  rrule::RefinementRule,
  space::FESpace,
  reffes::AbstractVector{<:ReferenceFE}
)
  cface_to_own_fface_to_own_dofs, _ = 
    _compute_cface_to_own_fface_to_own_dofs_and_permutations(rrule,space,reffes)
  return cface_to_own_fface_to_own_dofs
end

# This function computes at the same time the data structures
#    [coarse face][owned fine face] -> [owned dofs]
#    [coarse face][owned fine face][local permutation] -> [local permuted dofs]
# which then are used to compute the owned dofs and the permutations for each cface.
function _compute_cface_to_own_fface_to_own_dofs_and_permutations(
  rrule::RefinementRule{<:Polytope{Dc}},
  space::FESpace,
  reffes::AbstractVector{<:ReferenceFE};
) where Dc
  poly = get_polytope(rrule)
  topo = get_grid_topology(rrule.ref_grid)
  
  coffsets = get_offsets(poly)
  foffsets = get_offsets(topo)

  # Get RefinementRule topological information.
  # We use low-level functions to avoid repeating work
  d_to_fface_to_cface, d_to_fface_to_cface_dim = get_d_to_face_to_parent_face(rrule)
  cface_to_num_own_ffaces = _compute_cface_to_num_own_ffaces(
    rrule, d_to_fface_to_cface, d_to_fface_to_cface_dim
  )
  cface_to_own_ffaces = _compute_cface_to_own_ffaces(
    rrule, d_to_fface_to_cface, d_to_fface_to_cface_dim, cface_to_num_own_ffaces
  )
  d_to_cell_to_lface = map(Df -> Geometry.get_faces(topo,Dc,Df),0:Dc)

  cell_to_dofs = get_cell_dof_ids(space)
  cell_to_lface_to_dof = lazy_map(get_face_own_dofs,reffes)
  cell_to_lface_to_fpindex_to_ldofs = lazy_map(get_face_own_dofs_permutations,reffes)
  cell_to_offsets = lazy_map(r -> get_offsets(get_polytope(r)),reffes)

  # We need to collect the dofs for each fine face, but the owned dofs are given cell-wise.
  # So we need to iterate over the cells and identify which fine face we are looking at...
  # Since we are going to see some faces more than once, we keep track of which ones we 
  # have already seen through the `touched` array.
  # We also collect local dof permutations on that face, which will be aggregated later.
  cface_to_own_fface_to_dofs = [
    Vector{Vector{Int}}(undef,cface_to_num_own_ffaces[cface]) for cface in 1:num_faces(poly)
  ]
  cface_to_own_fface_to_fpindex_to_ldofs = [
    Vector{Vector{Vector{Int}}}(undef,cface_to_num_own_ffaces[cface]) for cface in 1:num_faces(poly)
  ]
  touched = fill(false,num_faces(topo))

  # For each subcell in the macro-element
  for cell in 1:num_subcells(rrule)
    dofs = view(cell_to_dofs,cell)
    lface_to_dof = cell_to_lface_to_dof[cell]
    lface_to_fpindex_to_ldofs = cell_to_lface_to_fpindex_to_ldofs[cell]
    offsets = cell_to_offsets[cell]
    for d in 0:Dc
      o = offsets[d+1] # Offset for the local d-faces
      fface_to_cface = d_to_fface_to_cface[d+1]
      fface_to_cface_dim = d_to_fface_to_cface_dim[d+1]
      # For each d-face of the subcell
      for (lface,fdface) in enumerate(d_to_cell_to_lface[d+1][cell])
        # fdface is the id of the fface (within the fine d-faces)
        fface = foffsets[d+1] + fdface # Id of the fine face (within the grid)
        if !touched[fface]
          face_dofs = view(dofs,lface_to_dof[o+lface])
          face_fpindex_to_ldofs = lface_to_fpindex_to_ldofs[o+lface]
          cdface = fface_to_cface[fdface]           # Id of the cface (within the coarse d-faces)
          cdface_dim = fface_to_cface_dim[fdface]   # Dimension of the cface
          cface = coffsets[cdface_dim+1] + cdface   # Id of the cface (within the polytope)
          pos = searchsortedfirst(cface_to_own_ffaces[cface],fface) # Id of fface (within the cface)
          cface_to_own_fface_to_dofs[cface][pos] = collect(Int,face_dofs)
          cface_to_own_fface_to_fpindex_to_ldofs[cface][pos] = face_fpindex_to_ldofs
          touched[fface] = true
        end
      end
    end
  end
  @check all(touched)

  return cface_to_own_fface_to_dofs, cface_to_own_fface_to_fpindex_to_ldofs
end

"""
    get_cface_to_own_dof_permutations(
      rrule::RefinementRule{<:Polytope},
      space::FESpace,
      reffes::AbstractVector{<:ReferenceFE}
    )

Given a RefinementRule and information on the dofs owned by each fine face, compute 
the permutations of the dofs owned by each coarse face.
"""
function get_cface_to_own_dof_permutations(
  rrule::RefinementRule{<:Polytope},
  space::FESpace,
  reffes::AbstractVector{<:ReferenceFE}
)
  # Dof permutation data
  cface_to_fface_to_dofs, 
    cface_to_fface_to_fpindex_to_ldofs = 
      _compute_cface_to_own_fface_to_own_dofs_and_permutations(rrule,space,reffes)

  # We need to convert the dof numberings (which are in the rrule numeration) into 
  # the local numeration of the coarse faces. 
  cface_to_fface_to_ldofs = map(cface_to_fface_to_dofs) do fface_to_dofs
    fface_to_ldofs = Vector{Vector{Int}}(undef,length(fface_to_dofs))
    offset = 0
    for (fface,dofs) in enumerate(fface_to_dofs)
      n_ldofs = length(dofs)
      fface_to_ldofs[fface] = collect(offset+1:offset+n_ldofs)
      offset += n_ldofs
    end
    return fface_to_ldofs
  end

  # Topological permutation data
  cface_to_cpindex_to_ffaces, 
    cface_to_cpindex_to_fpindex = get_cface_to_own_fface_permutations(rrule)

  # Allocate output
  poly = get_polytope(rrule)
  n_cfaces = num_faces(poly)
  n_cperms = map(length,cface_to_cpindex_to_ffaces)
  cface_to_cpindex_to_dofs = [
    Vector{Vector{Int}}(undef,n_cperms[cface]) for cface in 1:n_cfaces
  ]

  # For each coarse face, and each coarse permutation of the cface
  for cface in 1:n_cfaces
    fface_to_dofs = cface_to_fface_to_ldofs[cface]
    fface_to_fpindex_to_ldofs = cface_to_fface_to_fpindex_to_ldofs[cface]
    for cpindex in 1:n_cperms[cface]
      p_ffaces = cface_to_cpindex_to_ffaces[cface][cpindex]      # Permuted fine faces
      fface_pindex = cface_to_cpindex_to_fpindex[cface][cpindex] # Local fine permutation index

      # Collect dofs for each fine face owned by the coarse face, 
      # and permute them according to the two-level permutation induced by 
      # the coarse permutation (see `get_cface_to_own_fface_permutations`)
      dofs = Int32[]
      for (p_fface,fpindex) in zip(p_ffaces,fface_pindex)
        face_dofs = fface_to_dofs[p_fface]
        local_dof_permutation = fface_to_fpindex_to_ldofs[p_fface][fpindex]
        dofs = vcat(dofs,face_dofs[local_dof_permutation])
      end
      cface_to_cpindex_to_dofs[cface][cpindex] = dofs
    end
  end

  return cface_to_cpindex_to_dofs
end