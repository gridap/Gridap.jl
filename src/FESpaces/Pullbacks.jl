function get_cell_shapefuns_and_dof_basis(
  model::DiscreteModel, cell_reffe::AbstractArray{T}, conf::Conformity;
  scale_dof=false, global_meshsize=nothing
) where T <: ReferenceFE

  pushforward = Pushforward(get_name(T), conf)
  piola_map, cell_changes, cell_push_args = get_cell_pushforward(
    pushforward, model, cell_reffe)

  cell_ref_fields = lazy_map(get_shapefuns, cell_reffe)
  cell_ref_dofs = lazy_map(get_dof_basis, cell_reffe)

  if piola_map isa IdentityPiolaMap
    # The precomposition with geomap is handled by change_domain
    cell_phy_fields = cell_ref_fields
    cell_phy_dofs = cell_ref_dofs
  else
    @assert !isnothing(cell_changes) # needed to get a chance of applying DOFScalingMap next
    cell_phy_fields = lazy_map(piola_map, cell_ref_fields, cell_push_args...)
    dof_pf = inverse_map(Pullback(piola_map))
    cell_phy_dofs = lazy_map(dof_pf, cell_ref_dofs, cell_push_args...)
  end

  if isnothing(cell_changes)
    cell_phy_fields, cell_phy_dofs
  else
    cell_change, cell_change_invt = apply_dof_scaling(pushforward, model, cell_reffe, cell_changes..., scale_dof, global_meshsize)
    cell_shapefuns = lazy_map(linear_combination, cell_change, cell_phy_fields)
    cell_dof_basis = lazy_map(linear_combination, cell_change_invt, cell_phy_dofs)
    cell_shapefuns, cell_dof_basis
  end
end

function get_cell_pushforward(
  ::Pushforward, model::DiscreteModel, cell_reffe
)
  @abstractmethod
end

function get_cell_pushforward(
  p::IdentityPiolaMap, model::DiscreteModel, cell_reffe
)
  p, nothing, ()
end

# ContraVariantPiolaMap

function get_cell_pushforward(
  p::ContraVariantPiolaMap, model::DiscreteModel, cell_reffe
)
  cell_map  = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = get_sign_flip(model, cell_reffe) # equal to its transposed inverse
  return p, (change,change), (Jt,)
end

# CoVariantPiolaMap

function get_cell_pushforward(
  p::CoVariantPiolaMap, model::DiscreteModel, cell_reffe
)
  cell_map = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = lazy_map(r -> Diagonal(ones(num_dofs(r))), cell_reffe) # TODO: Replace by edge-signs
  return p, (change,change), (Jt,)
end

# DoubleContraVariantPiolaMap

using Gridap.ReferenceFEs: DoubleContraVariantPiolaMap
function get_cell_pushforward(
  ::DoubleContraVariantPiolaMap, model::DiscreteModel, cell_reffe, conformity
)
  cell_map = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = lazy_map(r -> Diagonal(fill(one(Float64), num_dofs(r))), cell_reffe) # TODO: Replace by sign flip
  #change  = get_sign_flip(model, cell_reffe) # equal to its transposed inverse
  return DoubleContraVariantPiolaMap(), (change,change), (Jt,)
end


#################
# NormalSignMap #
#################

function get_sign_flip(model::DiscreteModel{Dc}, cell_reffe) where Dc
  # Comment: lazy_maps on cell_reffes are very optimised, since they are CompressedArray/FillArray
  get_facet_own_dofs(reffe) = view(get_face_own_dofs(reffe),get_dimrange(get_polytope(reffe),Dc-1))
  cell_facet_own_dofs = lazy_map(get_facet_own_dofs, cell_reffe)
  cell_ids = IdentityVector(Int32(num_cells(model)))
  return lazy_map(NormalSignMap(model), cell_reffe, cell_facet_own_dofs, cell_ids)
end

"""
    struct NormalSignMap <: Map
      ...
    end

The `NormalSignMap` compute the signs to apply to the mapped reference normals,
for each facet of a physical cell.

Each physical facet ``f`` is shared by up to two cells ``K`` and ``K'``. The
orientation of the physical/global normal to ``f`` is chosen by the main cell
``K``, the first one in the list of adjascent cells to ``f`` in the grid topology.

The physical/global normal is the (normalized) Piola mapped reference normal to
``f̂ = F⁻¹(f)`` where ``F`` is the geometrical map ``F:K̂->K``. It is also minus
the (normalized) Piola mapped reference normal to ``f̂' = F'⁻¹(f)`` where ``F'``
is the geometrical map ``F':K̂->K'``.
"""
struct NormalSignMap{T} <: Map
  model::T
  facet_owners::Vector{Int32}
end

function NormalSignMap(model)
  facet_owners = compute_facet_owners(model)
  NormalSignMap(model,facet_owners)
end

function return_value(k::NormalSignMap,reffe,facet_own_dofs,cell)
  Diagonal(fill(one(Float64), num_dofs(reffe)))
end

function return_cache(k::NormalSignMap,reffe,facet_own_dofs,cell)
  model = k.model
  Dc = num_cell_dims(model)
  topo = get_grid_topology(model)

  cell_facets = get_faces(topo, Dc, Dc-1)
  cell_facets_cache = array_cache(cell_facets)

  return cell_facets, cell_facets_cache, CachedVector(Float64)
end

function evaluate!(cache,k::NormalSignMap,reffe,facet_own_dofs,cell)
  cell_facets,cell_facets_cache,dof_sign_cache = cache
  facet_owners = k.facet_owners

  setsize!(dof_sign_cache, (num_dofs(reffe),))
  dof_sign = dof_sign_cache.array

  o = one(eltype(dof_sign))
  fill!(dof_sign, o)

  facets = getindex!(cell_facets_cache,cell_facets,cell)
  for (lfacet,facet) in enumerate(facets)
    owner = facet_owners[facet]
    if owner != cell
      for dof in facet_own_dofs[lfacet]
        dof_sign[dof] = -o
      end
    end
  end

  return Diagonal(dof_sign)
end

function compute_facet_owners(model::DiscreteModel{Dc}) where {Dc}
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


#################
# DOFScalingMap #
#################

function apply_dof_scaling(pushforward, model, cell_reffe, change, change_invt,
                           scale_dof, global_meshsize)
  if scale_dof
    cell_ids = IdentityVector(Int32(num_cells(model)))
    dofscaling_map = DOFScalingMap(model, cell_reffe, pushforward, global_meshsize)
    dof_scales = lazy_map(dofscaling_map, cell_ids)
    change = lazy_map(*, change, dof_scales)
    inv_dof_scales = lazy_map(inv, dof_scales)
    change_invt = lazy_map(*, inv_dof_scales, change_invt)
  end
  change, change_invt
end

"""
    struct DOFScalingMap{M,S,V} <: Map
"""
struct DOFScalingMap{M,S,V} <: Map
  d_id_to_data::M
  scaler::S
  cell_ctype::V
  ctype_ndofs::Vector{Int}
  ctype_offsets::Vector{Vector{Int}}
  ctype_faceowndofs::Vector{Vector{Vector{Int}}}

  @doc """
      DOFScalingMap(model, cell_reffe, pushforward, global_meshsize=nothing)


  """
  function DOFScalingMap(model, cell_reffe, pushforward, ::Nothing)
    ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
    ctype_ndof = num_dofs.(ctype_reffe)
    ctype_offsets = @. get_offsets(get_polytope(ctype_reffe))
    ctype_faceowndofs = get_face_own_dofs.(ctype_reffe)

    fdims = unique_dim_of_faces_owning_dofs(ctype_reffe)

    Dc = num_cell_dims(model)
    topo = get_grid_topology(model)
    d_id_to_data = ()
    for d in fdims
      cell_to_dfaces = get_faces(topo, Dc, d)
      cell_to_dfaces_cache = array_cache(cell_to_dfaces)
      dface_to_fmeas = get_cell_measure(Triangulation(ReferenceFE{d},model))
      dface_to_fmeas_cache = array_cache(dface_to_fmeas)
      d_data = (d, cell_to_dfaces, cell_to_dfaces_cache, dface_to_fmeas, dface_to_fmeas_cache)
      d_id_to_data = (d_id_to_data..., d_data)
    end

    scaler = dof_scaling_function(pushforward, Dc)

    new{typeof(d_id_to_data), typeof(scaler), typeof(cell_ctype)}(
      d_id_to_data, scaler, cell_ctype, ctype_ndof, ctype_offsets, ctype_faceowndofs
    )
  end

  # version with prescribed global_meshsize
  # the d_id_to_data field is re-purposed to hold global_meshsize
  function DOFScalingMap(model, cell_reffe, pushforward, global_meshsize)
    @assert global_meshsize isa Number
    ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
    ctype_ndof = num_dofs.(ctype_reffe)
    Dc = num_cell_dims(model)
    scaler = dof_scaling_function(pushforward, Dc)

    new{typeof(global_meshsize), typeof(scaler), typeof(cell_ctype)}(
      global_meshsize, scaler, cell_ctype, ctype_ndof, [Int[]], [[Int[]]]
    )
  end
end

function return_value(s::DOFScalingMap, cell)
  cell_type = s.cell_ctype[cell]
  @inbounds ndofs = s.ctype_ndofs[cell_type]
  Diagonal(fill(one(Float64), ndofs))
end

function return_cache(s::DOFScalingMap, cell)
  cell_type = s.cell_ctype[cell]
  @inbounds ndofs = s.ctype_ndofs[cell_type]
  cache = CachedVector(Float64)
  setsize!(cache, (ndofs,))
  cache
end

function evaluate!(cache, s::DOFScalingMap, cell)
  cell_type = s.cell_ctype[cell]
  @inbounds ndofs = s.ctype_ndofs[cell_type]
  @inbounds offsets = s.ctype_offsets[cell_type]
  @inbounds face_own_dofs = s.ctype_faceowndofs[cell_type]

  setsize!(cache, (ndofs,))
  dof_scale = cache.array
  @check begin fill!(dof_scale, 0); true end # to check that all scales are set later

  function _set_d_faces_dofscales!(dof_scale, cell, d_data)
      d, c_to_dfaces, c_to_dfaces_cache, f_to_fmeas, f_to_fmeas_cache = d_data
      @inbounds dfaces = getindex!(c_to_dfaces_cache, c_to_dfaces, cell)
      @inbounds dfaces_offset = offsets[d+1]

      @inbounds for (lface, face) in enumerate(dfaces)
        rface = dfaces_offset + lface # face index in polytope
        isempty(face_own_dofs[rface]) && continue
        face_dvol = getindex!(f_to_fmeas_cache, f_to_fmeas, face)
        h = face_dvol^(1/d)
        face_dofscale = s.scaler(h)
        for dof in face_own_dofs[rface]
          dof_scale[dof] = face_dofscale
        end
      end
  end

  map(data -> _set_d_faces_dofscales!(dof_scale, cell, data), s.d_id_to_data)

  @check all(!iszero, dof_scale) "Some DOF scale have not been set. Open an issue, and/or disable DOF scaling"
  return Diagonal(dof_scale)
end

# version with prescribed global_meshsize
function evaluate!(cache, s::DOFScalingMap{<:Number}, cell)
  cell_type = s.cell_ctype[cell]
  @inbounds ndofs = s.ctype_ndofs[cell_type]

  setsize!(cache, (ndofs,))
  dof_scale = cache.array

  global_meshsize = s.d_id_to_data
  dofscale = s.scaler(global_meshsize)
  @check !iszero(dofscale)
  dof_scale .= dofscale

  return Diagonal(dof_scale)
end


"""
    unique_dim_of_faces_owning_dofs(ctype_reffe) -> Vector{Int}

Vector of all unique dimension of a face owning a DOF in reffes in `ctype_reffe`.
"""
function unique_dim_of_faces_owning_dofs(ctype_reffe)
  dims_owning_dof = Int[]
  for reffe in ctype_reffe
    face_own_dofs = get_face_own_dofs(reffe)
    owning_faces = findall(!isempty, face_own_dofs)
    dimranges = get_dimranges(get_polytope(reffe))
    reffe_dims = map(f -> findfirst(∋(f), dimranges) - 1, owning_faces) |> unique!
    union!(dims_owning_dof, reffe_dims)
  end
  dims_owning_dof
end

