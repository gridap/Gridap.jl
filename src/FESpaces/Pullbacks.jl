
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
# Edge orientation signs     #
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

function compute_cell_bases_changes(
  ::Argyris, ::IdentityPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
)
  p = get_polytope(testitem(cell_reffe))
  cell_σ = _edge_signs(model, p)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  cell_change = lazy_map(
    ReferenceFEs.ArgyrisChangeOfBasis(p, false), cell_Jtx, cell_σ
  )
  cell_change_invt = lazy_map(
    ReferenceFEs.ArgyrisChangeOfBasis(p, true), cell_Jtx, cell_σ
  )
  return (cell_change, cell_change_invt)
end

############################################################################################
# Morley

function compute_cell_bases_changes(
  ::Morley, ::IdentityPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
)
  p = get_polytope(testitem(cell_reffe))
  cell_σ = _edge_signs(model, p)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  cell_change = lazy_map(
    ReferenceFEs.MorleyChangeOfBasis(p, false), cell_Jtx, cell_σ
  )
  cell_change_invt = lazy_map(
    ReferenceFEs.MorleyChangeOfBasis(p, true), cell_Jtx, cell_σ
  )
  return (cell_change, cell_change_invt)
end

############################################################################################
# HHJ, Regge and GLS

# These three share `EdgeScalingChangeOfBasis` verbatim: each has an edge DoF of
# the form ∫ₑ (a⋅Mb) μᵢ ds for a pair of directions carried dually by that
# element's push-forward, so the Jacobians cancel and only 1/‖J t̂ₑ‖ is left. See
# that type, in src/ReferenceFEs/Pullbacks.jl, for why the three coincide. The
# interior DoFs are left as the push-forward of the reference ones -- they are
# cell-owned and shared with nobody, so any per-cell convention gives the same
# space, as Gridap already does for the cell moments of Raviart-Thomas and BDM.

function _edge_scaling_cell_bases_changes(model::DiscreteModel, cell_reffe, cell_Jt)
  reffe = testitem(cell_reffe)
  p = get_polytope(reffe)
  cell_σ = _edge_signs(model, p)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  cell_change = lazy_map(
    ReferenceFEs.EdgeScalingChangeOfBasis(reffe, false), cell_Jtx, cell_σ
  )
  cell_change_invt = lazy_map(
    ReferenceFEs.EdgeScalingChangeOfBasis(reffe, true), cell_Jtx, cell_σ
  )
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

function compute_cell_bases_changes(
  ::ArnoldWintherNC, ::ReferenceFEs.DoubleContraVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt
)
  K = ReferenceFEs.AWNCChangeOfBasis
  _aw_cell_bases_changes(K, model, cell_reffe, cell_Jt)
end

function compute_cell_bases_changes(
  ::ArnoldWintherC, ::ReferenceFEs.DoubleContraVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt
)
  K = ReferenceFEs.AWCChangeOfBasis
  _aw_cell_bases_changes(K, model, cell_reffe, cell_Jt)
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

function compute_cell_bases_changes(
  ::MardalTaiWinther, ::ContraVariantPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
)
  reffe = testitem(cell_reffe)
  p = get_polytope(reffe)

  # The geometrical map is affine on simplices, so its Jacobian is constant.
  x0 = Fill(first(get_vertex_coordinates(p)), length(cell_Jt))
  cell_Jtx = lazy_map(evaluate, cell_Jt, x0)

  if num_dims(p) == 2
    cell_σ = _edge_signs(model, p)
    cell_change = lazy_map(
      ReferenceFEs.MTWChangeOfBasis(p, false), cell_Jtx, cell_σ
    )
    cell_change_invt = lazy_map(
      ReferenceFEs.MTWChangeOfBasis(p, true), cell_Jtx, cell_σ
    )
  else
    cell_pids = get_cell_permutations(get_grid_topology(model), 2)
    cell_change = lazy_map(
      ReferenceFEs.TWChangeOfBasis(reffe, false), cell_Jtx, cell_pids
    )
    cell_change_invt = lazy_map(
      ReferenceFEs.TWChangeOfBasis(reffe, true), cell_Jtx, cell_pids
    )
  end
  return (cell_change, cell_change_invt)
end

############################################################################################
# Cartesian products

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
  return (lazy_map(ReferenceFEs.CartProdBlockDiag(), fwd...),
          lazy_map(ReferenceFEs.CartProdBlockDiag(), bwd...))
end

function compute_cell_bases_changes(
  ::CartProd{K,NS}, pf::ReferenceFEs.CartProdPushforward{K},
  model::DiscreteModel, cell_reffe, cell_Jt
) where {K,NS}
  pfs = ReferenceFEs._cp_maps(pf, Val(K))
  _cp_cell_bases_changes(Val(K), pfs, model, cell_reffe, cell_Jt)
end

# every factor is identity-mapped, so the stack is too and the push-forward comes
# back unwrapped
function compute_cell_bases_changes(
  ::CartProd{K,NS}, pf::IdentityPiolaMap, model::DiscreteModel, cell_reffe, cell_Jt
) where {K,NS}
  pfs = ReferenceFEs._cp_maps(pf, Val(K))
  _cp_cell_bases_changes(Val(K), pfs, model, cell_reffe, cell_Jt)
end
