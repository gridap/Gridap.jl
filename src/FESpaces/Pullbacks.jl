
function get_cell_shapefuns_and_dof_basis(
  model::DiscreteModel, cell_reffe::AbstractArray{T}, conf::Conformity; kwargs...
) where T <: ReferenceFE

  cell_map  = get_cell_map(get_grid(model))
  cell_Jt = lazy_map(Broadcasting(∇), cell_map)

  reffe_name = get_name(T)
  pushforward = Pushforward(reffe_name, conf)
  cell_changes = compute_cell_bases_changes(reffe_name, pushforward, model, cell_reffe, cell_Jt)

  return get_cell_shapefuns_and_dof_basis(
    pushforward, model, cell_reffe, cell_changes, cell_Jt; kwargs...
  )
end

# This constructor allows fo provided cell changes and jacobians, 
# which is necessary for GridapDistributed
function get_cell_shapefuns_and_dof_basis(
  pushforward::Pushforward, model, cell_reffe, cell_changes, cell_Jt;
  scale_dof=false, global_meshsize=nothing
)
  cell_ref_fields = lazy_map(get_shapefuns, cell_reffe)
  cell_ref_dofs = lazy_map(get_dof_basis, cell_reffe)

  # Apply the pushforward "individually" to each shape-function, and the inverse pullback to each DOF
  if pushforward isa IdentityPiolaMap
    # The precomposition with geomap is handled by change_domain
    cell_phy_fields = cell_ref_fields
    cell_phy_dofs = cell_ref_dofs
  else
    cell_phy_fields = lazy_map(pushforward, cell_ref_fields, cell_Jt)
    dof_pf = inverse_map(Pullback(pushforward))
    cell_phy_dofs = lazy_map(dof_pf, cell_ref_dofs, cell_Jt)
  end

  # If nontrivial, apply the appropriate change of basis to the DOF and shape-function bases
  cell_changes = apply_dof_scaling(
    cell_changes, model, cell_reffe, pushforward, scale_dof, global_meshsize
  )
  if isnothing(cell_changes)
    return (cell_phy_fields, cell_phy_dofs)
  end

  cell_change, cell_change_invt = cell_changes
  cell_shapefuns = lazy_map(linear_combination, cell_change,      cell_phy_fields)
  cell_dof_basis = lazy_map(linear_combination, cell_change_invt, cell_phy_dofs)
  return cell_shapefuns, cell_dof_basis
end

"""
    compute_cell_bases_changes(name::ReferenceFEName, push::Pushforward, model, cell_reffe, cell_Jt)

Computes, in each cell, the change of basis ``M`` between the pushforwarded
reference shape-function basis and the expected physical shape-functions, as
well as it transposed inverse ``M⁻ᵀ``, the change of basis between the inverse
pullback of the reference cell DOF and the expected physical cell DOF.

For that, the `model` is provided. As well as the reference FE + lazy gradient
(transposed Jacobian) of the geometrical map, in each cell.

See also the manual page on ["FE basis transformations"](@ref "FE basis transformations").
Return either `nothing` (no change required) or a couple of cell arrays `(cell_M, cell_M⁻ᵀ)`

The `dof_scale` and `global_meshsize` kwargs are not handeled by this function.
For them to work properly, it might be necessary to also overload
[`get_dofscale_setter_function`](@ref).
"""
function compute_cell_bases_changes(
  name::ReferenceFEName, push::Pushforward, model, cell_reffe, cell_Jt
)
  @abstractmethod
end

function compute_cell_bases_changes(
  ::ReferenceFEName, ::IdentityPiolaMap, model, cell_reffe, cell_Jt
)
  nothing
end

function compute_cell_bases_changes(
  ::ReferenceFEName, ::ContraVariantPiolaMap, model, cell_reffe, cell_Jt
)
  change = get_sign_flip(model, cell_reffe) # equal to its transposed inverse
  return (change,change)
end

function compute_cell_bases_changes(
  ::ReferenceFEName, ::CoVariantPiolaMap, model, cell_reffe, cell_Jt
)
  D = num_cell_dims(model)
  poly = only(get_polytopes(model))
  if (D==2) || is_simplex(poly)
    # For these cases, we do not need to aply a sign flip
    return nothing
  elseif (D==3) && is_n_cube(poly)
    change = get_sign_flip(model, cell_reffe)
    return (change,change)
  end
  @notimplemented
end

using Gridap.ReferenceFEs: DoubleContraVariantPiolaMap
function compute_cell_bases_changes(
  ::ReferenceFEName, ::DoubleContraVariantPiolaMap, model, cell_reffe, cell_Jt
)
  change = lazy_map(r -> Diagonal(fill(one(Float64), num_dofs(r))), cell_reffe) # TODO: Replace by sign flip
  #change  = get_sign_flip(model, cell_reffe) # equal to its transposed inverse
  return (change,change)
end

#################
# NormalSignMap #
#################

function get_sign_flip(model::DiscreteModel, cell_reffe, sign_map = NormalSignMap(model))
  # Comment: lazy_maps on cell_reffes are very optimised, since they are CompressedArray/FillArray
  Dc = num_cell_dims(model)
  get_facet_own_dofs(reffe) = view(get_face_own_dofs(reffe),get_dimrange(get_polytope(reffe),Dc-1))
  cell_facet_own_dofs = lazy_map(get_facet_own_dofs, cell_reffe)
  cell_ids = IdentityVector(Int32(num_cells(model)))
  return lazy_map(sign_map, cell_reffe, cell_facet_own_dofs, cell_ids)
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

function compute_facet_owners(model::DiscreteModel{Dc}, select_nbor=maximum) where {Dc}
  topo = get_grid_topology(model)
  facet_to_cell = get_faces(topo, Dc-1, Dc)

  nfacets = num_faces(topo, Dc-1)
  owners = Vector{Int32}(undef, nfacets)
  for facet in 1:nfacets
    facet_cells = view(facet_to_cell, facet)
    @check !isempty(facet_cells) "Facet $facet has no adjacent cells"
    selected_owner = select_nbor(facet_cells)
    @check selected_owner isa Integer "select_nbor must return an integer owner for facet $facet, got $(typeof(selected_owner))"
    owner = Int(selected_owner)
    @check owner != 0 "select_nbor returned invalid owner 0 for facet $facet; expected one of $(collect(facet_cells))"
    @check owner in facet_cells "select_nbor returned invalid owner $owner for facet $facet; expected one of $(collect(facet_cells))"
    owners[facet] = Int32(owner)
  end
  return owners
end

##############################
# EdgeScalingChangeOfBasis   #
##############################

function _edge_signs(model::DiscreteModel, p::Polytope{2})
  cell_ledge_pindex = get_cell_permutations(get_grid_topology(model), 1)
  nedges = num_faces(p, 1)

  cache = array_cache(cell_ledge_pindex)
  signs = Vector{NTuple{nedges,Float64}}(undef, length(cell_ledge_pindex))
  for cell in eachindex(cell_ledge_pindex)
    pinds = getindex!(cache, cell_ledge_pindex, cell)
    signs[cell] = ntuple(e -> ifelse(isone(pinds[e]), 1.0, -1.0), nedges)
  end
  return signs
end

# The change of basis of an element whose edge DoFs are *preserved up to one
# positive scalar per edge*: diagonal, with `‖J t̂ₑ‖^{±1}` on the DoFs of edge `e`,
# times `(-1)ⁱ` when the cell traverses that edge against the global direction, and
# `1` on the interior DoFs. `transposed_inverse` selects `P⁻ᵀ` over `P`.
#
# This covers every element whose edge DoF has the form `∫ₑ (a⋅Mb) μᵢ ds` for a
# pair of directions carried dually by that element's push-forward — a normal and a
# normal under the double contravariant map, two tangents under the double
# covariant one, one of each under the co-contravariant one. In every case the
# Jacobians cancel completely,
#
#     a⋅Mb = (â⋅M̂b̂) / ‖J t̂ₑ‖²,   ds = ‖J t̂ₑ‖ dŝ   ⟹   F∗(ℓ^{e,i}) = ℓ̂^{e,i} / ‖J t̂ₑ‖.
#
# That is not a coincidence. A Piola map is *chosen* so that its element's DoF is
# invariant, so what is left over cannot depend on which pairing was picked — only
# on the fact that the DoF is an edge moment against a degree-`i` weight. The
# elements still need their own `compute_cell_bases_changes`, which dispatches on
# the reference FE name and the push-forward.
#
# Since `a⋅Mb` is quadratic in the directions, or bilinear with both of them
# flipping, none of these elements needs a normal or tangent sign convention; the
# only orientation effect is the parity of the weight, which is the `(-1)ⁱ` above.
# At order 0 even that disappears.
struct EdgeScalingChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  edge_dofs::Vector{Vector{Int}}
  ndofs::Int
  transposed_inverse::Bool
end

function EdgeScalingChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  edge_dofs = [own[nv+e] for e in 1:num_faces(p, 1)]
  return EdgeScalingChangeOfBasis(
    get_edge_tangent(p), edge_dofs, num_dofs(reffe), transposed_inverse
  )
end

function return_cache(k::EdgeScalingChangeOfBasis, Jt, σ)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs))
end

function evaluate!(cache, k::EdgeScalingChangeOfBasis, Jt, σ)
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  for i in 1:k.ndofs
    M[i, i] = 1.0
  end

  for e in eachindex(k.edge_dofs)
    L = norm(k.tangents[e] ⋅ Jt)   # ‖J t̂ₑ‖, the edge length ratio
    reversed = σ[e] < 0
    for (i, d) in enumerate(k.edge_dofs[e])
      # the i-th DoF of the edge carries the weight of degree i-1
      s = ifelse(reversed && isodd(i - 1), -1.0, 1.0)
      M[d, d] = ifelse(k.transposed_inverse, s / L, s * L)
    end
  end

  return M
end

#################
# DOFScalingMap #
#################

function apply_dof_scaling(cell_changes, model, cell_reffe, pushforward,
                           scale_dof, global_meshsize)

  !scale_dof  && return cell_changes
  # scaling Dof is necessary if either:
  # - the piola map is non-IdentityPiolaMap,
  # - or the cell change is non-trivial (e.g. C1 reffes, that will still use IdentityPiolaMap)
  isnothing(cell_changes) && pushforward isa IdentityPiolaMap && return cell_changes

  if isnothing(cell_changes)
    cell_change = lazy_map(r -> Diagonal(fill(one(Float64), num_dofs(r))), cell_reffe)
    cell_change_invt = cell_change
  else
    cell_change, cell_change_invt = cell_changes
  end

  cell_ids = IdentityVector(Int32(num_cells(model)))
  dofscaling_map = DOFScalingMap(model, cell_reffe, pushforward, global_meshsize)
  dof_scales = lazy_map(dofscaling_map, cell_ids)
  inv_dof_scales = lazy_map(inv, dof_scales)

  cell_change = lazy_map(*, dof_scales, cell_change)
  cell_change_invt = lazy_map(*, inv_dof_scales, cell_change_invt)
  return (cell_change, cell_change_invt)
end

struct DOFScalingMap{M,V,S} <: Map
  # this is global_meshsize::Real if given,
  # else data on cell faces and faces volumes for each dimension d of faces owning DOFs in any reffe.
  d_id_to_data::M

  # data related to the unique reffes ≡ ctype in the given cell_reffe
  cell_ctype::V
  ctype_scalesetter::S # function(s) returned by get_dofscale_setter_function
  ctype_ndofs::Vector{Int}
  ctype_offsets::Vector{Vector{Int}}
  ctype_faceowndofs::Vector{Vector{Vector{Int}}}

  @doc """
      DOFScalingMap(model, cell_reffe, pushforward, nothing)
      DOFScalingMap(model, cell_reffe, pushforward, global_meshsize::Real)

  Evaluated at a `cell` id of `model` to return a `Diagonal` matrix that rescales
  the cell physical shape-function basis, which is the result of mapping the
  `cell_reffe[cell]` basis using `pushforward`.

  `global_meshsize` is either `nothing` or a `Real` number. If `nothing`, the
  local meshsize is estimated on each `d` dimensional face using the `d`-root
  of its `d`-volume.

  `DOFScalingMap` is designed to be robust to heterogeneous cell reffes, and to
  reffes having heterogeneous DOF scaling (even within a face) like
  Mardal-Tai-Winter or C1 reffes.

  The method to overload in order to implement nontrivial `DOFScalingMap` to a
  new reffe is [`get_dofscale_setter_function`](@ref).
  See also ["FE basis transformations"](@ref "FE basis transformations").
  """
  function DOFScalingMap(model, cell_reffe, pushforward, ::Nothing)
    ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
    ctype_ndof = num_dofs.(ctype_reffe)
    ctype_offsets = @. get_offsets(get_polytope(ctype_reffe))
    ctype_faceowndofs = get_face_own_dofs.(ctype_reffe)
    ctype_scalesetter = Tuple(
      get_dofscale_setter_function(reffe, pushforward) for reffe in ctype_reffe
    )

    fdims = unique_dim_of_faces_owning_dofs(ctype_reffe)
    @notimplementedif 0 in fdims """
    No local meshsize estimator at physical node is currently implemented, but a reffe has vertices owned DOF.
    It is currently only possible to use `scale_dof` with `global_meshsize` in this case.

    The reason is that the meshsize is currently estimated as the `d`-root of the `d`-volume, which is always 1 for a vertex.
    """

    Dc = num_cell_dims(model)
    topo = get_grid_topology(model)
    d_id_to_data = ()
    for d in fdims
      if iszero(d)
        # Options:
        # - use cell diameter or D-root of cell volume,
        # of a master cell picked like in NormalSignMap, or averaged over all adjascent cells
      else
        cell_to_dfaces = get_faces(topo, Dc, d)
        cell_to_dfaces_cache = array_cache(cell_to_dfaces)
        dface_to_fmeas = get_cell_measure(Triangulation(ReferenceFE{d},model))
        dface_to_fmeas_cache = array_cache(dface_to_fmeas)
        d_data = (d, cell_to_dfaces, cell_to_dfaces_cache, dface_to_fmeas, dface_to_fmeas_cache)
      end
      d_id_to_data = (d_id_to_data..., d_data)
    end

    new{typeof(d_id_to_data), typeof(cell_ctype), typeof(ctype_scalesetter)}(
      d_id_to_data, cell_ctype, ctype_scalesetter, ctype_ndof, ctype_offsets, ctype_faceowndofs
    )
  end

  # version with prescribed global_meshsize
  # the d_id_to_data field is re-purposed to hold global_meshsize
  function DOFScalingMap(model, cell_reffe, pushforward, global_meshsize::Real)
    @assert global_meshsize isa Number
    ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
    ctype_ndof = num_dofs.(ctype_reffe)
    ctype_faceowndofs = get_face_own_dofs.(ctype_reffe)
    ctype_scalesetter = Tuple(
      get_dofscale_setter_function(reffe, pushforward) for reffe in ctype_reffe
    )

    new{typeof(global_meshsize), typeof(cell_ctype), typeof(ctype_scalesetter)}(
      global_meshsize, cell_ctype, ctype_scalesetter, ctype_ndof, [Int[]], ctype_faceowndofs
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
  @inbounds face_own_dofs = s.ctype_faceowndofs[cell_type]

  dof_scale_cache = CachedVector(Float64)
  setsize!(dof_scale_cache , (ndofs,))
  face_meshsize_cache = CachedVector(Float64)
  setsize!(face_meshsize_cache , (length(face_own_dofs),))

  dof_scale_cache, face_meshsize_cache
end

function evaluate!(cache, s::DOFScalingMap, cell)
  cell_type = s.cell_ctype[cell]
  @inbounds ndofs = s.ctype_ndofs[cell_type]
  @inbounds offsets = s.ctype_offsets[cell_type]
  @inbounds face_own_dofs = s.ctype_faceowndofs[cell_type]
  @inbounds scale_setter! = s.ctype_scalesetter[cell_type]

  dof_scale_cache, face_meshsize_cache = cache
  setsize!(dof_scale_cache, (ndofs,))
  dof_scale = dof_scale_cache.array
  setsize!(face_meshsize_cache, (length(face_own_dofs),))
  face_meshsize = face_meshsize_cache.array
  @check begin fill!(dof_scale, 0); true end # to check that all scales are set later
  @check begin fill!(face_meshsize, 0); true end

  function _set_d_faces_meshsizes!(face_meshsize, cell, d_data)
      d, c_to_dfaces, c_to_dfaces_cache, f_to_fmeas, f_to_fmeas_cache = d_data
      if iszero(d)
        @notimplemented # TODO
      else
        @inbounds dfaces = getindex!(c_to_dfaces_cache, c_to_dfaces, cell)
        @inbounds dfaces_offset = offsets[d+1]

        @inbounds for (lface, face) in enumerate(dfaces)
          rface = dfaces_offset + lface # face index in polytope
          isempty(face_own_dofs[rface]) && continue
          face_dvol = getindex!(f_to_fmeas_cache, f_to_fmeas, face)
          h = face_dvol^(1/d)
          face_meshsize[rface] = h
        end
    end
  end

  map(data -> _set_d_faces_meshsizes!(face_meshsize, cell, data), s.d_id_to_data)
  scale_setter!(dof_scale, face_own_dofs, face_meshsize)

  @check all(!iszero, dof_scale) "Some DOF scale have not been set. Open an issue, and/or disable DOF scaling"
  return Diagonal(dof_scale)
end

# version with prescribed global_meshsize
function evaluate!(cache, s::DOFScalingMap{<:Real}, cell)
  cell_type = s.cell_ctype[cell]
  @inbounds ndofs = s.ctype_ndofs[cell_type]
  @inbounds face_own_dofs = s.ctype_faceowndofs[cell_type]
  @inbounds scale_setter! = s.ctype_scalesetter[cell_type]

  dof_scale_cache, face_meshsize_cache = cache
  setsize!(dof_scale_cache, (ndofs,))
  dof_scale = dof_scale_cache.array
  setsize!(face_meshsize_cache, (length(face_own_dofs),))
  face_meshsize = face_meshsize_cache.array
  @check begin fill!(dof_scale, 0); true end # to check that all scales are set later

  global_meshsize = s.d_id_to_data
  fill!(face_meshsize, global_meshsize)
  scale_setter!(dof_scale, face_own_dofs, face_meshsize)

  @check all(!iszero, dof_scale) "Some DOF scale have not been set. Open an issue, and/or disable DOF scaling"
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

############################################################################################
# Rotating PΛ: per-cell change of basis
#
# Mesh-level conformity layer for the rotating P_rΛ¹ and trimmed P_r⁻Λ¹ reference
# FEs (see src/ReferenceFEs/RotatingPLambdaRefFEs.jl).  Each cell is re-expressed
# in the frame of its "virtually sorted" copy, at the permutation
#
#   π_K = sortperm(cell global vertex ids)
#
# so that two cells sharing a face agree on that face's dof order.  The matrices
# themselves come from `compute_pλ_change` (Gridap.Polynomials), which is pure
# rotation calculus and knows nothing about the mesh; the derivation of the
# convention lives there.
#
# Since π_K is a permutation of the D+1 cell vertices, at most (D+1)! distinct
# changes of basis can occur on any mesh — 6 in 2D, 24 in 3D, independently of
# the number of cells.  They are therefore built once, on demand, and shared
# through a `CompressedArray` indexed by a per-cell permutation id.
#
# Fast path: when every cell already lists its vertices in increasing global-id
# order (e.g. simplexified Cartesian models) π_K is the identity everywhere and
# no change of basis is needed at all.

function compute_cell_bases_changes(::RotatingPΛName, ::CoVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt)
  compute_pλ_cell_bases_changes(model, cell_reffe)
end

function compute_cell_bases_changes(::TrimmedPΛName, ::CoVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt)
  compute_pλ_cell_bases_changes(model, cell_reffe)
end

function compute_pλ_cell_bases_changes(model::DiscreteModel, cell_reffe)
  D    = num_cell_dims(model)
  topo = get_grid_topology(model)
  cell_verts = Geometry.get_faces(topo, D, 0)

  # One reference basis shared by all cells (single-reffe meshes).
  basis = get_prebasis(testitem(cell_reffe))
  rc    = RotationCache(basis)
  return compute_pλ_cell_bases_changes(cell_verts, rc)
end

function compute_pλ_cell_bases_changes(cell_verts, rc::RotationCache)
  cell_to_pid = zeros(Int8, length(cell_verts))
  pid_to_M = Matrix{Float64}[]
  pid_to_Minv = Matrix{Float64}[]

  alltrivial = true
  π_to_pid = Dict{Vector{Int},Int8}()
  cache = array_cache(cell_verts)
  for cell in eachindex(cell_verts)
    v = getindex!(cache, cell_verts, cell)
    π = sortperm(v)
    pid = get!(π_to_pid, π) do
      M, Minv = compute_pλ_change(rc, π)
      push!(pid_to_M, M)
      push!(pid_to_Minv, Minv)
      return Int8(length(pid_to_M))
    end
    alltrivial &= issorted(v)
    cell_to_pid[cell] = pid
  end

  # Fast path: every cell already lists its vertices in increasing global-id
  # order (e.g. simplexified Cartesian models) — π_K = id everywhere.
  alltrivial && return nothing

  # General case
  cell_change = CompressedArray(pid_to_M, cell_to_pid)
  cell_change_invt = CompressedArray(pid_to_Minv, cell_to_pid)
  return (cell_change, cell_change_invt)
end

############################################################################################
# Argyris

#     _congruence_matrix(A) -> TensorValue{3,3}
# 
# The 3×3 matrix of the congruence `H ↦ A H Aᵀ` acting on symmetric 2×2 matrices,
# in the coordinates `(H₁₁, H₁₂, H₂₂)` — equivalently the second symmetric power
# `Sym²(A)`, i.e. `A ⊗ A` restricted to the symmetric subspace.
function _congruence_matrix(A)
  a11, a12, a21, a22 = A[1,1], A[1,2], A[2,1], A[2,2]
  # column-major, one column per basis matrix: [1 0;0 0], [0 1;1 0], [0 0;0 1]
  return TensorValue{3,3}(
    a11*a11, a11*a21,             a21*a21,
    2*a11*a12, a11*a22 + a12*a21, 2*a21*a22,
    a12*a12, a12*a22,             a22*a22
  )
end

# Argyris is mapped by the plain pullback u = û∘F⁻¹. Writing K = J⁻ᵀ, and using
# that F is affine on a simplex so no second derivative of F appears,
#
#   ∇u = K ∇û,        D²u = K D²û Kᵀ,
#
# the vertex value DoFs are preserved, the vertex gradient DoFs mix within a
# vertex through K, and the vertex Hessian DoFs mix within a vertex through
# `_congruence_matrix(K)`. The edge DoFs behave as in Morley: with G = (JᵀJ)⁻¹
# and n = R t,
#
#   F∗(δₑᵏ) = Aₖ δ̂ₑᵏ + Bₖ (δ̂ᵥᵇ - δ̂ᵥᵃ),   Aₖ = det(J) n̂ᵀGn̂,  Bₖ = det(J) t̂ᵀGn̂,
#
# coupling each edge row to the two *value* DoFs of its endpoints. So
#
#         ⎡ I    0    0    0 ⎤                  ⎡ I       0     0      0   ⎤
#   W  =  ⎢ 0    Kg   0    0 ⎥ ,      W⁻¹  =    ⎢ 0       Kg⁻¹  0      0   ⎥
#         ⎢ 0    0    Kh   0 ⎥                  ⎢ 0       0     Kh⁻¹   0   ⎥
#         ⎣ Bv   0    0    A ⎦                  ⎣ -A⁻¹Bv  0     0      A⁻¹ ⎦
#
# with Kg, Kh block diagonal over the vertices and A = diag(Aₖ) — block
# triangular, not block diagonal. Kg⁻¹ and Kh⁻¹ are the same constructions
# applied to K⁻¹ = Jᵀ, so no matrix is ever inverted numerically.
#
# Edge orientation follows `_edge_signs`, with D = diag(1,…,1,σ₁,σ₂,σ₃) folded in
# as P = W⁻¹D and P⁻ᵀ = WᵀD. The vertex DoFs need no such treatment, being stated
# in the global Cartesian frame and so the same functional for every cell that
# touches the vertex.

#     ArgyrisChangeOfBasis(p, transposed_inverse)
#
# Builds, from the (transposed) Jacobian of a cell's geometrical map and that
# cell's edge orientation signs `σ`, either the change of basis `P`
# (`transposed_inverse = false`) or `P⁻ᵀ` (`transposed_inverse = true`).
struct ArgyrisChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  edge_vertices::Vector{Vector{Int}}
  transposed_inverse::Bool
end

function ArgyrisChangeOfBasis(p::Polytope{2}, transposed_inverse::Bool)
  ts, ns = ReferenceFEs._edge_frames(p)
  ArgyrisChangeOfBasis(ts, ns, get_faces(p, 1, 0), transposed_inverse)
end

function return_cache(k::ArgyrisChangeOfBasis, Jt, σ)
  nv, ne = length(k.tangents), length(k.edge_vertices)
  CachedArray(zeros(Float64, 6*nv + ne, 6*nv + ne))
end

function evaluate!(cache, k::ArgyrisChangeOfBasis, Jt, σ)
  nv = length(k.tangents)   # a triangle: as many vertices as edges
  ne = length(k.edge_vertices)
  ndofs = 6*nv + ne
  setsize!(cache, (ndofs, ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))

  detJ = det(Jt)
  G = inv(Jt ⋅ transpose(Jt))   # (JᵀJ)⁻¹, since Jt = Jᵀ

  # W carries K = J⁻ᵀ on the gradients and `_congruence_matrix(K)` on the
  # Hessians, so W⁻¹ carries K⁻¹ = Jᵀ = Jt and `_congruence_matrix(Jt)`, while Wᵀ
  # carries their transposes.
  K = inv(Jt)
  Kg = ifelse(k.transposed_inverse, transpose(K), Jt)
  Kh = k.transposed_inverse ? transpose(_congruence_matrix(K)) :
       _congruence_matrix(Jt)

  gof = nv           # offset of the gradient DoFs
  hof = 3*nv         # offset of the Hessian DoFs
  eof = 6*nv         # offset of the edge DoFs

  for v in 1:nv
    M[v, v] = 1.0
    for i in 1:2, j in 1:2
      M[gof+2*v-2+i, gof+2*v-2+j] = Kg[i, j]
    end
    for i in 1:3, j in 1:3
      M[hof+3*v-3+i, hof+3*v-3+j] = Kh[i, j]
    end
  end

  for e in 1:ne
    t̂, n̂ = k.tangents[e], k.normals[e]
    Gn̂ = G ⋅ n̂
    A = detJ * (n̂ ⋅ Gn̂)
    B = detJ * (t̂ ⋅ Gn̂)
    σe = σ[e]
    va, vb = k.edge_vertices[e]

    # D = diag(1,…,1,σ…) scales the last columns of the block below.
    if k.transposed_inverse
      # WᵀD
      M[eof+e, eof+e] = A * σe
      M[va, eof+e] = -B * σe
      M[vb, eof+e] = B * σe
    else
      # W⁻¹D
      M[eof+e, eof+e] = σe / A
      M[eof+e, va] = B / A
      M[eof+e, vb] = -B / A
    end
  end

  return M
end

function compute_cell_bases_changes(
  ::Argyris, ::IdentityPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
)
  p = get_polytope(testitem(cell_reffe))
  cell_σ = _edge_signs(model, p)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  cell_change = lazy_map(ArgyrisChangeOfBasis(p, false), cell_Jtx, cell_σ)
  cell_change_invt = lazy_map(ArgyrisChangeOfBasis(p, true), cell_Jtx, cell_σ)
  return (cell_change, cell_change_invt)
end

############################################################################################
# Morley

# Morley is mapped by the plain pullback u = û∘F⁻¹, under which the vertex values
# are preserved but the edge normal derivatives are not: with ∇u = J⁻ᵀ∇û,
#
#   ∇u⋅n = ∇û⋅(J⁻¹n) = a (∇û⋅n̂) + b (∇û⋅t̂),
#
# so the push-forward of an edge DoF picks up a tangential derivative, which is
# not a Morley node. With the DoFs written as moments that tangential piece is
# exactly a difference of vertex values,
#
#   ∫_{ê} ∇û⋅t̂ dŝ = û(v̂_b) - û(v̂_a),
#
# so span(N̂) is preserved after all and the transformation is a 6×6 matrix in
# closed form. With n = R t and G = (JᵀJ)⁻¹, using R J Rᵀ = det(J) J⁻ᵀ,
#
#   J⁻¹n = (det J / ‖J t̂‖) G n̂,   ds = ‖J t̂‖ dŝ,
#
# the ‖J t̂‖ cancels and, for edge k with reference endpoints v̂_a, v̂_b,
#
#   F∗(δᵥⁱ) = δ̂ᵥⁱ,
#   F∗(δₑᵏ) = Aₖ δ̂ₑᵏ + Bₖ (δ̂ᵥᵇ - δ̂ᵥᵃ),   Aₖ = det(J) n̂ᵀGn̂,  Bₖ = det(J) t̂ᵀGn̂.
#
# So W = [I 0; B A] with A = diag(Aₖ) — block *triangular*, the case Kirby notes
# for Morley and Argyris. Aₖ ≠ 0 always (n̂ᵀGn̂ > 0 by positive definiteness), and
# W⁻¹ = [I 0; -A⁻¹B A⁻¹] in closed form. Edge orientation follows `_edge_signs`,
# with D = diag(1,1,1,σ₁,σ₂,σ₃) folded in as P = W⁻¹D and P⁻ᵀ = WᵀD.

#     MorleyChangeOfBasis(p, transposed_inverse)
#
# Builds, from the (transposed) Jacobian of a cell's geometrical map and that
# cell's edge orientation signs `σ`, either the change of basis `P`
# (`transposed_inverse = false`) or `P⁻ᵀ` (`transposed_inverse = true`).
struct MorleyChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  edge_vertices::Vector{Vector{Int}}
  transposed_inverse::Bool
end

function MorleyChangeOfBasis(p::Polytope{2}, transposed_inverse::Bool)
  ts, ns = ReferenceFEs._edge_frames(p)
  MorleyChangeOfBasis(ts, ns, get_faces(p, 1, 0), transposed_inverse)
end

function return_cache(k::MorleyChangeOfBasis, Jt, σ)
  ndofs = length(k.tangents) + length(k.edge_vertices)
  CachedArray(zeros(Float64, ndofs, ndofs))
end

function evaluate!(cache, k::MorleyChangeOfBasis, Jt, σ)
  nedges = length(k.tangents)
  nverts = nedges  # a triangle
  ndofs = nverts + nedges
  setsize!(cache, (ndofs, ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))

  detJ = det(Jt)
  G = inv(Jt ⋅ transpose(Jt))  # (JᵀJ)⁻¹, since Jt = Jᵀ

  for i in 1:nverts
    M[i, i] = 1.0
  end

  for e in 1:nedges
    t̂, n̂ = k.tangents[e], k.normals[e]
    Gn̂ = G ⋅ n̂
    A = detJ * (n̂ ⋅ Gn̂)
    B = detJ * (t̂ ⋅ Gn̂)
    σe = σ[e]
    va, vb = k.edge_vertices[e]

    # D = diag(1,1,1,σ...) scales the last columns of the block below.
    if k.transposed_inverse
      # WᵀD
      M[nverts+e, nverts+e] = A * σe
      M[va, nverts+e] = -B * σe
      M[vb, nverts+e] = B * σe
    else
      # W⁻¹D = [I 0; -A⁻¹B A⁻¹D]
      M[nverts+e, nverts+e] = σe / A
      M[nverts+e, va] = B / A
      M[nverts+e, vb] = -B / A
    end
  end

  return M
end

function compute_cell_bases_changes(
  ::Morley, ::IdentityPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
)
  p = get_polytope(testitem(cell_reffe))
  cell_σ = _edge_signs(model, p)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  cell_change = lazy_map(MorleyChangeOfBasis(p, false), cell_Jtx, cell_σ)
  cell_change_invt = lazy_map(MorleyChangeOfBasis(p, true), cell_Jtx, cell_σ)
  return (cell_change, cell_change_invt)
end

############################################################################################
# HHJ, Regge and GLS

# These three share `EdgeScalingChangeOfBasis` verbatim: each has an edge DoF of
# the form ∫ₑ (a⋅Mb) μᵢ ds for a pair of directions carried dually by that
# element's push-forward, so the Jacobians cancel and only 1/‖J t̂ₑ‖ is left. See
# that type above for why the three coincide. The interior DoFs are left as the
# push-forward of the reference ones -- they are cell-owned and shared with
# nobody, so any per-cell convention gives the same space, as Gridap already does
# for the cell moments of Raviart-Thomas and BDM.

function _edge_scaling_cell_bases_changes(model::DiscreteModel, cell_reffe, cell_Jt)
  reffe = testitem(cell_reffe)
  p = get_polytope(reffe)
  cell_σ = _edge_signs(model, p)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  cell_change = lazy_map(EdgeScalingChangeOfBasis(reffe, false), cell_Jtx, cell_σ)
  cell_change_invt = lazy_map(EdgeScalingChangeOfBasis(reffe, true), cell_Jtx, cell_σ)
  return (cell_change, cell_change_invt)
end

function compute_cell_bases_changes(
  ::HellanHerrmannJohnson, ::ReferenceFEs.DoubleContraVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt
)
  _edge_scaling_cell_bases_changes(model, cell_reffe, cell_Jt)
end

function compute_cell_bases_changes(
  ::Regge, ::ReferenceFEs.DoubleCoVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt
)
  _edge_scaling_cell_bases_changes(model, cell_reffe, cell_Jt)
end

function compute_cell_bases_changes(
  ::GopalakrishnanLedererSchoberl, ::ReferenceFEs.CoContraVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt
)
  _edge_scaling_cell_bases_changes(model, cell_reffe, cell_Jt)
end

############################################################################################
# Arnold-Winther, conforming and nonconforming

# Both AW elements are mapped by the double contravariant Piola map
# τ = det(J)⁻² J τ̂ Jᵀ, and both carry the same four edge DoFs -- the degree 0 and
# 1 moments of n⋅τn and of n⋅τt. Those give one 4×4 block per edge, in the DoF
# order (nn0, nn1, nt0, nt1),
#
#   W = [1/L  0    0  0;  0  1/L  0  0;  α  0  β  0;  0  α  0  β],
#
# with L = ‖J t̂ₑ‖, α = ã/(det(J) L), β = L/det(J) and ã = n̂ᵀ(JᵀJ)t̂. The interior
# DoFs are left as the push-forward of the reference ones -- they are cell-owned
# and shared with nobody -- which sidesteps a dense interior block entirely, as
# Gridap already does for the cell moments of Raviart-Thomas.
#
# The conforming element adds a 3×3 block per vertex, det(J)⁻² times the matrix
# of H ↦ J H Jᵀ, i.e. `_congruence_matrix` with A = J rather than the A = J⁻ᵀ the
# Argyris Hessian block uses. Both DoF kinds are invariant under reversing an
# edge -- n⋅τn is quadratic in n, n⋅τt bilinear with both flipping -- so only the
# parity of the Legendre weight enters σ.

function _aw_edge_blocks!(M, tangents, normals, edge_dofs, Jt, σ, transposed_inverse)
  detJ = det(Jt)
  JtJ = Jt ⋅ transpose(Jt)   # JᵀJ, since Jt = Jᵀ

  for e in eachindex(edge_dofs)
    t̂, n̂ = tangents[e], normals[e]
    L = norm(t̂ ⋅ Jt)               # ‖J t̂‖
    α = (n̂ ⋅ (JtJ ⋅ t̂)) / (detJ * L)
    β = L / detJ
    reversed = σ[e] < 0

    # DoF order within an edge: (nn,0), (nn,1), (nt,0), (nt,1)
    dofs = edge_dofs[e]
    nmom = length(dofs) ÷ 2
    for i in 1:nmom
      s = ifelse(reversed && isodd(i - 1), -1.0, 1.0)   # parity of the Legendre weight
      dnn, dnt = dofs[i], dofs[nmom+i]
      if transposed_inverse
        # WᵀD
        M[dnn, dnn] = s / L
        M[dnn, dnt] = s * α
        M[dnt, dnt] = s * β
      else
        # W⁻¹D
        M[dnn, dnn] = s * L
        M[dnt, dnn] = -s * α * L / β
        M[dnt, dnt] = s / β
      end
    end
  end
end

#     AWNCChangeOfBasis(reffe, transposed_inverse)
#
# One 4×4 block per edge, and the identity on the interior.
struct AWNCChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  edge_dofs::Vector{Vector{Int}}
  ndofs::Int
  transposed_inverse::Bool
end

function AWNCChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  ts, ns = ReferenceFEs._edge_frames(p)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  edge_dofs = [own[nv+e] for e in 1:num_faces(p, 1)]
  AWNCChangeOfBasis(ts, ns, edge_dofs, num_dofs(reffe), transposed_inverse)
end

function return_cache(k::AWNCChangeOfBasis, Jt, σ)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs))
end

function evaluate!(cache, k::AWNCChangeOfBasis, Jt, σ)
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  for i in 1:k.ndofs
    M[i, i] = 1.0     # the interior DoFs are left as the push-forward
  end
  _aw_edge_blocks!(M, k.tangents, k.normals, k.edge_dofs, Jt, σ, k.transposed_inverse)
  return M
end

#     AWCChangeOfBasis(reffe, transposed_inverse)
#
# A 3×3 block per vertex, a 4×4 block per edge, and the identity on the interior.
struct AWCChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  vertex_dofs::Vector{Vector{Int}}
  edge_dofs::Vector{Vector{Int}}
  ndofs::Int
  transposed_inverse::Bool
end

function AWCChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  ts, ns = ReferenceFEs._edge_frames(p)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  vertex_dofs = [own[v] for v in 1:nv]
  edge_dofs = [own[nv+e] for e in 1:num_faces(p, 1)]
  AWCChangeOfBasis(ts, ns, vertex_dofs, edge_dofs, num_dofs(reffe), transposed_inverse)
end

function return_cache(k::AWCChangeOfBasis, Jt, σ)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs))
end

function evaluate!(cache, k::AWCChangeOfBasis, Jt, σ)
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  for i in 1:k.ndofs
    M[i, i] = 1.0     # the interior DoFs are left as the push-forward
  end

  J = transpose(Jt)
  detJ = det(Jt)

  # W carries det(J)⁻² _congruence_matrix(J) on the vertices, so W⁻¹ carries
  # det(J)² _congruence_matrix(J⁻¹) and Wᵀ the transpose of the former.
  Bv = k.transposed_inverse ? transpose(_congruence_matrix(J)) / detJ^2 :
       _congruence_matrix(inv(J)) * detJ^2
  for dofs in k.vertex_dofs
    for i in 1:3, j in 1:3
      M[dofs[i], dofs[j]] = Bv[i, j]
    end
  end

  _aw_edge_blocks!(M, k.tangents, k.normals, k.edge_dofs, Jt, σ, k.transposed_inverse)
  return M
end

function compute_cell_bases_changes(
  ::ArnoldWintherNC, ::ReferenceFEs.DoubleContraVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt
)
  _aw_cell_bases_changes(AWNCChangeOfBasis, model, cell_reffe, cell_Jt)
end

function compute_cell_bases_changes(
  ::ArnoldWintherC, ::ReferenceFEs.DoubleContraVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt
)
  _aw_cell_bases_changes(AWCChangeOfBasis, model, cell_reffe, cell_Jt)
end

function _aw_cell_bases_changes(K, model::DiscreteModel, cell_reffe, cell_Jt)
  reffe = testitem(cell_reffe)
  p = get_polytope(reffe)
  cell_σ = _edge_signs(model, p)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  cell_change = lazy_map(K(reffe, false), cell_Jtx, cell_σ)
  cell_change_invt = lazy_map(K(reffe, true), cell_Jtx, cell_σ)
  return (cell_change, cell_change_invt)
end

############################################################################################
# Mardal-Tai-Winther, 2D and 3D
#
# The element is one family but its change of basis has two shapes, because the
# mismatch between two cells sharing a facet is a *sign* in 2D and a
# *permutation* in 3D: an edge has two orderings, a triangular face has 3! = 6.
# `compute_cell_bases_changes` therefore branches on the polytope, feeding
# `_edge_signs` in 2D and `get_cell_permutations(topo, 2)` in 3D. Both blocks are
# closed form; neither inverts anything at run time.

# 2D: block diagonal, one 3x3 block per edge, acting on (ℓⁿ⁰, ℓⁿ¹, ℓᵗ⁰). Normal
# moments are preserved by the contravariant Piola map, while tangential ones
# pick up a normal component:
#
#   Wᵏ = [s 0 0; 0 s 0; α 0 β],  s = sign(det J),
#   α = (J n̂)⋅(J t̂)/|det J|,     β = ‖J t̂‖²/|det J|,
#
# post-multiplied by Dᵏ = diag(σᵏ,1,σᵏ). This is
# [Aznaran, Farrell & Kirby, SMAI-JCM 8 (2022) 399, (5.21)-(5.25)]; the
# degree-1 normal moment contributes an identity row, the parity of μ₁ being what
# keeps ℓⁿ¹ from mixing into ℓⁿ⁰.
struct MTWChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  normals::Vector{VectorValue{2,Float64}}
  transposed_inverse::Bool
end

function return_cache(k::MTWChangeOfBasis, Jt, σ)
  ndofs = 3*length(k.tangents)
  CachedArray(zeros(Float64, ndofs, ndofs))
end

function evaluate!(cache, k::MTWChangeOfBasis, Jt, σ)
  nedges = length(k.tangents)
  setsize!(cache, (3*nedges, 3*nedges))
  M = cache.array
  fill!(M, zero(eltype(M)))

  s = sign(det(Jt))
  adetJ = abs(det(Jt))

  for e in 1:nedges
    # Jt is the transposed Jacobian, so `v ⋅ Jt` is `J v`.
    Jt̂ = k.tangents[e] ⋅ Jt
    Jn̂ = k.normals[e] ⋅ Jt
    α = (Jn̂ ⋅ Jt̂) / adetJ
    β = (Jt̂ ⋅ Jt̂) / adetJ
    σe = σ[e]

    # Dᵏ = diag(σ,1,σ) scales the columns of the block below.
    o = 3*(e - 1)
    if k.transposed_inverse
      # (Wᵏ)ᵀDᵏ
      M[o+1, o+1] = s * σe
      M[o+2, o+2] = s
      M[o+1, o+3] = α * σe
      M[o+3, o+3] = β * σe
    else
      # (Wᵏ)⁻¹Dᵏ
      M[o+1, o+1] = s * σe
      M[o+2, o+2] = s
      M[o+3, o+1] = -s * α * σe / β
      M[o+3, o+3] = σe / β
    end
  end

  return M
end

# 3D: one 6x6 block per face, and nothing else.
#
# The block is closed form. Its DoF weights are *equivariant* -- every one is a
# fixed reference object pushed forward by the map, `n = e₁×e₂` with `eᵢ = J êᵢ`,
# and the tangential weights `n×e₁`, `n×e₂`, `n×(x-a)` -- so the pulled-back
# physical weights expand in the reference weights with constant coefficients:
#
#   Jᵀwⱼ = [detJ·I₃  0 ;  Ã  λI₃] ŵ,   λ = |n|²/|n̂|²,  a = (Jn̂)⋅(n×e₁)/|n̂|²
#
# (and `b` likewise with `e₂`), whence `W = that / detJ`. The normal moments come
# out preserved exactly -- the identity block -- and `W⁻¹` is closed form too, so
# no quadrature, no cached shape-function moments and no matrix inversion appear
# at run time.
#
# This is the 3D analogue of [AFK, (5.21)-(5.25)], which they give only in 2D;
# FIAT/FInAT define the 3D element but transform it only on triangles. Getting it
# needs the equivariant weight basis: with an intrinsic frame (normalised
# tangents, unit normal) the coefficients stop being constant and the block no
# longer closes without integrating the shape functions.
struct TWChangeOfBasis <: Map
  face_dofs::Vector{Vector{Int}}          # dofs owned by each face
  # per (face, permutation), the reference edge vectors of the face taken in the
  # mesh's own vertex order -- the only reference data the block needs
  ref_edges::Matrix{NTuple{2,VectorValue{3,Float64}}}
  # per (face, permutation), `inv(T)` or `transpose(T)` as needed, where `T`
  # re-expresses the reference weights of the mesh-sorted frame in those of the
  # polytope's own frame -- see `_tw_frame_change`
  pre::Matrix{Matrix{Float64}}
  ndofs::Int
  transposed_inverse::Bool
end

function TWChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  own = get_face_own_dofs(reffe)
  fdim = get_dimrange(p, 2)
  face_dofs = [own[f] for f in fdim]
  nfaces = length(face_dofs)

  fverts = get_faces(p, 2, 0)
  vcoords = get_vertex_coordinates(p)
  vperms = get_face_vertex_permutations(p, 2)   # the list Gridap's pindex indexes

  # the same P⁻₁Λ¹ (RT₀) weight basis the tangential moments are declared with
  μb = FEEC_poly_basis(Val(2), Float64, 1, 1, :P⁻; rotate_90=true)
  nouts = get_facet_normal(p)                   # the outward normals `σ` uses
  nperms = length(first(vperms))
  ref_edges = Matrix{NTuple{2,VectorValue{3,Float64}}}(undef, nfaces, nperms)
  pre = Matrix{Matrix{Float64}}(undef, nfaces, nperms)
  for f in 1:nfaces, (pid, perm) in enumerate(vperms[f])
    # `pindex` maps cell-face-vertex -> global-face-vertex, so its inverse lists
    # the cell's local vertices in the order the mesh stores the face
    lv = fverts[f][invperm(perm)]
    â, b̂, ĉ = vcoords[lv[1]], vcoords[lv[2]], vcoords[lv[3]]
    ref_edges[f, pid] = (b̂ - â, ĉ - â)

    loc = fverts[f]
    T = _tw_frame_change((vcoords[loc[1]], vcoords[loc[2]], vcoords[loc[3]]),
                         (â, b̂, ĉ), μb, perm, nouts[f])
    pre[f, pid] = transposed_inverse ? transpose(T) : inv(T)
  end

  TWChangeOfBasis(face_dofs, ref_edges, pre, num_dofs(reffe), transposed_inverse)
end

# The 6x6 matrix `T` with `ŵˢⱼ = Σₖ Tⱼₖ ŵˡₖ`, expressing the canonical weights of
# the mesh-sorted frame of a face in the ones the moments actually declare, in
# the polytope's own frame. Both blocks are closed form.
#
# The normal block is `ε A`, with `ε = ±1` the sign the scaled normal `n̂ₛ = ê₁×ê₂`
# picks up under the relabelling and `A` the affine change from `(1, uˢ, vˢ)` to
# `(1, uˡ, vˡ)`, determined exactly by three points.
#
# The tangential block is `Mμ⁻¹ R`, where `R = rotation_change_of_basis(μb,
# invperm(π))` is the vertex relabelling of `μ` and `Mμ` holds the `RT₀`
# coefficients `(a, b, c)` of `μⱼ = (a + c u, b + c v)`, which is the change from
# `{Ĵf μⱼ}` to the canonical `{ê₁, ê₂, x̂-â}`. `μ` is a trimmed FEEC space, so `R`
# is available in closed form and, at this degree, is a signed permutation with
# integer entries. See Badia, Manyer & Marteau, *Rotating bases for finite
# element exterior calculus*.
#
# Both factors are exact, so nothing here is fitted. `Mμ` is geometry-free and
# `R` is integer; the only reason `T` is assembled per (face, permutation) at all
# is the normal block's dependence on the face's own vertex labelling.
function _tw_frame_change(lverts, sverts, μb, π, nout)
  v3(w) = Float64[w[1], w[2], w[3]]
  T = zeros(6, 6)

  a, e1, e2 = lverts[1], lverts[2] - lverts[1], lverts[3] - lverts[1]
  as, es1, es2 = sverts[1], sverts[2] - sverts[1], sverts[3] - sverts[1]
  ns = cross(e1, e2)
  ε = ns ⋅ cross(es1, es2) > 0 ? 1.0 : -1.0
  # the moments state their weights with the polytope's *outward* normal, which
  # carries no orientation of the shared face. `sout` is the sign relating it to
  # ê₁×ê₂, the orientation the sorted frame does carry, and `area2` undoes the
  # extension's normalisation of the tangential weights. Both are per-face
  # constants.
  sout = nout ⋅ ns > 0 ? 1.0 : -1.0
  area2 = norm(ns)
  Es = hcat(v3(es1), v3(es2))

  # (1, uˢ, vˢ) against (1, uˡ, vˡ) at three points fixes the affine relabelling
  Vl = zeros(3, 3)
  Vs = zeros(3, 3)
  for (i, (u, v)) in enumerate(((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)))
    uv = Es \ v3((a + u * e1 + v * e2) - as)
    Vl[i, :] = [1.0, u, v]
    Vs[i, :] = [1.0, uv[1], uv[2]]
  end
  T[1:3, 1:3] = sout * ε * transpose(Vl \ Vs)

  # μⱼ = (a + c u, b + c v), so Ĵf μⱼ = a ê₁ + b ê₂ + c (x̂ - â)
  m00 = evaluate(μb, [Point(0.0, 0.0)])
  m10 = evaluate(μb, [Point(1.0, 0.0)])
  Mμ = zeros(3, 3)
  for j in 1:3
    Mμ[j, :] = [m00[1, j][1], m00[1, j][2], m10[1, j][1] - m00[1, j][1]]
  end
  T[4:6, 4:6] = sout * area2 *
                inv(Mμ) * Polynomials.rotation_change_of_basis(μb, invperm(collect(π)))
  T
end

function return_cache(k::TWChangeOfBasis, Jt, pids)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs)), zeros(6, 6)
end

function evaluate!(_cache, k::TWChangeOfBasis, Jt, pids)
  cache, W = _cache
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  detJt = det(Jt)
  adetJ = abs(detJt)
  sgn = sign(detJt)

  for f in eachindex(k.face_dofs)
    ê1, ê2 = k.ref_edges[f, pids[f]]
    n̂ = cross(ê1, ê2)
    s2 = n̂ ⋅ n̂
    e1 = ê1 ⋅ Jt                        # `v ⋅ Jt` is `J v`
    e2 = ê2 ⋅ Jt
    n = cross(e1, e2)
    Jn̂ = n̂ ⋅ Jt

    μ = (n ⋅ n) / (s2 * adetJ)
    a = (Jn̂ ⋅ cross(n, e1)) / (s2 * adetJ)
    b = (Jn̂ ⋅ cross(n, e2)) / (s2 * adetJ)

    fill!(W, 0.0)
    if k.transposed_inverse                       # Wᵀ
      for i in 1:3
        W[i, i] = sgn
        W[3+i, 3+i] = μ
      end
      W[1, 4] = a; W[1, 5] = b; W[2, 6] = a; W[3, 6] = b
    else                                          # W⁻¹
      iμ = 1 / μ
      for i in 1:3
        W[i, i] = sgn
        W[3+i, 3+i] = iμ
      end
      W[4, 1] = -sgn * a * iμ; W[5, 1] = -sgn * b * iμ
      W[6, 2] = -sgn * a * iμ; W[6, 3] = -sgn * b * iμ
    end

    # the physical dofs are stated in the mesh-sorted frame, the reference ones
    # in the polytope's own, so the frame change composes in
    P = k.pre[f, pids[f]]
    dofs = k.face_dofs[f]
    for i in 1:6, j in 1:6
      M[dofs[i], dofs[j]] = sum(P[i, m] * W[m, j] for m in 1:6)
    end
  end

  return M
end

function compute_cell_bases_changes(
  ::MardalTaiWinther, ::ContraVariantPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
)
  reffe = testitem(cell_reffe)
  p = get_polytope(reffe)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  if num_dims(p) == 2
    ts, ns = ReferenceFEs._edge_frames(p)
    cell_σ = _edge_signs(model, p)
    cell_change = lazy_map(MTWChangeOfBasis(ts, ns, false), cell_Jtx, cell_σ)
    cell_change_invt = lazy_map(MTWChangeOfBasis(ts, ns, true), cell_Jtx, cell_σ)
  else
    cell_pids = get_cell_permutations(get_grid_topology(model), 2)
    cell_change = lazy_map(TWChangeOfBasis(reffe, false), cell_Jtx, cell_pids)
    cell_change_invt = lazy_map(TWChangeOfBasis(reffe, true), cell_Jtx, cell_pids)
  end
  return (cell_change, cell_change_invt)
end

############################################################################################
# Cartesian products

# The factors of a stacked element never mix, so its change of basis is the
# factors' assembled block diagonally, in the blocked ordering of
# `CartProdRefFEs.jl`. A factor that needs no change of basis contributes an
# identity block; when none of them do, `nothing` propagates.

struct CartProdBlockDiag <: Map end

function return_cache(::CartProdBlockDiag, Ps::AbstractMatrix...)
  T = promote_type(map(eltype, Ps)...)
  CachedArray(zeros(T, sum(P -> size(P, 1), Ps), sum(P -> size(P, 2), Ps)))
end

function evaluate!(cache, ::CartProdBlockDiag, Ps::AbstractMatrix...)
  setsize!(cache, (sum(P -> size(P, 1), Ps), sum(P -> size(P, 2), Ps)))
  A = cache.array
  fill!(A, zero(eltype(A)))
  io = jo = 0
  for P in Ps
    ni, mj = size(P)
    @inbounds A[io+1:io+ni, jo+1:jo+mj] .= P
    io += ni
    jo += mj
  end
  return A
end

_cp_maps(::ReferenceFEs.CartProdPushforward{K,PFS}, ::Val{K}) where {K,PFS} =
  ReferenceFEs._cp_maps(PFS)

_cp_maps(pf::IdentityPiolaMap, ::Val{K}) where K = ntuple(i -> pf, K)

function _cp_cell_bases_changes(::Val{K}, pfs, model, cell_reffe, cell_Jt) where K
  reffes = ReferenceFEs.get_metadata(testitem(cell_reffe))
  ncells = length(cell_reffe)

  changes = ntuple(K) do c
    compute_cell_bases_changes(
      get_name(reffes[c]), pfs[c], model, Fill(reffes[c], ncells), cell_Jt
    )
  end
  all(isnothing, changes) && return nothing

  eye(c) = Fill(Matrix{Float64}(I, num_dofs(reffes[c]), num_dofs(reffes[c])), ncells)
  fwd = ntuple(c -> isnothing(changes[c]) ? eye(c) : changes[c][1], K)
  bwd = ntuple(c -> isnothing(changes[c]) ? eye(c) : changes[c][2], K)
  return (lazy_map(CartProdBlockDiag(), fwd...),
          lazy_map(CartProdBlockDiag(), bwd...))
end

function compute_cell_bases_changes(
  ::CartProd{K,NS}, pf::ReferenceFEs.CartProdPushforward{K},
  model::DiscreteModel, cell_reffe, cell_Jt
) where {K,NS}
  _cp_cell_bases_changes(Val(K), _cp_maps(pf, Val(K)), model, cell_reffe, cell_Jt)
end

# every factor is identity-mapped, so the stack is too and the push-forward comes
# back unwrapped
function compute_cell_bases_changes(
  ::CartProd{K,NS}, pf::IdentityPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
) where {K,NS}
  _cp_cell_bases_changes(Val(K), _cp_maps(pf, Val(K)), model, cell_reffe, cell_Jt)
end
