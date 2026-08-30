# RotatingPLambda/PΛCellBasesChanges.jl
#
# Mesh-level conformity layer for the rotating P_rΛ¹ and trimmed P_r⁻Λ¹
# reference FEs (see src/ReferenceFEs/RotatingPLambdaRefFEs.jl): the
# compute_cell_bases_changes methods providing the per-cell change of basis
# given by the closed-form vertex-relabelling calculus of Gridap.Polynomials.
#
# ─────────────────────────────────────────────────────────────────────────────
# Per-cell change of basis: the rotation calculus at π_K = sortperm(vertex ids)
# ─────────────────────────────────────────────────────────────────────────────
#
# CONVENTION (derived below; pinned empirically over all relative vertex
# orderings of two-triangle and two-tet meshes — the competing candidates
# C(π)ᵀ, C(invperm(π)) and every un-relabelled variant fail the tangential
# jump/interpolation round-trip tests at O(1)):
#
# Let g be the cell's global vertex ids in local order and π = π_K =
# sortperm(g), so π maps sorted position → local index. The virtually sorted
# geometric map is F̃ = F ∘ A⁻¹ with A = A_{invperm(π)} (reference
# automorphism V_j ↦ V_{invperm(π)[j]}), hence the conforming basis is
#
#   φ̃_μ = (F̃⁻¹)^* ŵ_μ = (F⁻¹)^* (A_{invperm(π)}^* ŵ_μ)
#        = Σ_ν C(π)[μ,ν] · (F⁻¹)^* ŵ_ν,
#
# by the pullback identity
# (A_{invperm(π)}^* w_μ)(x) = Σ_ν C[μ,ν] w_ν(x) with
# C = rotation_change_of_basis(b, π).
#
# RELABELLING: φ̃_μ is supported on the LOCAL face π(F(μ)) (its bubble face
# key F(μ) is a virtual — rank-relabelled — key), while the conformity
# machinery glues local index μ through the reffe's face_own_dofs, i.e. by the
# local face F(μ). The local shape function stored at index μ must therefore
# be φ̃_{s(μ)}, where s is the face-wise order-preserving reindexing that maps
# the i-th bubble of local face F to the i-th bubble of the virtual face
# sort(invperm(π)(F)) — so that support(φ̃_{s(μ)}) = π(F(s(μ))) = F(μ), and
# every cell adjacent to a geometric face enumerates that face's functions in
# the bubble order of its sorted global vertex ids (this shared order is what
# makes the identity face-own-dof permutations correct for every pindex).
#
# The mathematical change matrix is M[μ,ν] = C(π)[s(μ),ν], and in the
# linear_combination convention (out[j] = Σ_i values[i,j]·in[i]):
#
#   cell_change[i,j]      = M[j,i]    = C(π)[s(j), i],
#   cell_change_invt[i,j] = M⁻¹[i,j]  = C(invperm(π))[i, s(j)]
#
# (using C(π)·C(invperm(π)) = I), the latter being the dof-side correction
# that restores duality. Both are assembled from the memoised sparse
# rotation_map rows — no matrix inversion.

using Gridap.Polynomials: RotationCache, rotation_map
using Gridap.Arrays: IdentityVector

struct _PΛRotationCoBMap{T,RC<:RotationCache} <: Map
  cell_verts  :: T
  rc          :: RC
  invert      :: Bool  # false → cell_change; true → cell_change_invt
  face_to_ws  :: Dict{Vector{Int},Vector{Int}}  # bubble face key → its w's (in order)
  memo        :: Dict{Vector{Int},Matrix{Float64}}
end

function _PΛRotationCoBMap(cell_verts, rc::RotationCache, invert::Bool)
  face_to_ws = Dict{Vector{Int},Vector{Int}}()
  for (w, (F, _, _)) in enumerate(rc.entries)
    push!(get!(face_to_ws, F, Int[]), w)
  end
  _PΛRotationCoBMap(cell_verts, rc, invert, face_to_ws,
                    Dict{Vector{Int},Matrix{Float64}}())
end

# s(μ): face-wise order-preserving reindexing local entry → virtual entry,
# mapping the i-th bubble of face F to the i-th bubble of sort(invperm(π)(F)).
function _pλ_sorted_relabel(face_to_ws, π::Vector{Int})
  invπ = invperm(π)
  n = sum(length(ws) for ws in values(face_to_ws))
  s = Vector{Int}(undef, n)
  for (F, ws) in face_to_ws
    wsG = face_to_ws[sort(invπ[F])]
    for i in eachindex(ws)
      s[ws[i]] = wsG[i]
    end
  end
  s
end

# NOTE: qualified extension — inside FESpaces the bare name `return_value`
# is a module-local function (created in Pullbacks.jl), not Arrays.return_value.
function Arrays.return_value(m::_PΛRotationCoBMap, cell)
  n = length(m.rc.entries)
  Matrix{Float64}(LinearAlgebra.I, n, n)
end

return_cache(m::_PΛRotationCoBMap, cell) = array_cache(m.cell_verts)

function evaluate!(cache, m::_PΛRotationCoBMap, cell)
  gverts = getindex!(cache, m.cell_verts, cell)
  π = sortperm(gverts)
  get!(m.memo, π) do
    n = length(m.rc.entries)
    s = _pλ_sorted_relabel(m.face_to_ws, π)
    M = zeros(Float64, n, n)
    if m.invert
      # M[w, s⁻¹(w′)] += C(invperm(π))[w, w′]
      sinv = invperm(s)
      rows = rotation_map(m.rc, invperm(π))
      for w in 1:n, (c, w′) in rows[w]
        M[w, sinv[w′]] += c
      end
    else
      # M[w′, μ] += C(π)[s(μ), w′]
      rows = rotation_map(m.rc, π)
      for μ in 1:n, (c, w′) in rows[s[μ]]
        M[w′, μ] += c
      end
    end
    M
  end
end

function compute_cell_bases_changes(::RotatingPΛName, ::CoVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt)
  _pλ_cell_bases_changes(model, cell_reffe)
end

function compute_cell_bases_changes(::TrimmedPΛName, ::CoVariantPiolaMap,
  model::DiscreteModel, cell_reffe, cell_Jt)
  _pλ_cell_bases_changes(model, cell_reffe)
end

function _pλ_cell_bases_changes(model::DiscreteModel, cell_reffe)
  D    = num_cell_dims(model)
  topo = get_grid_topology(model)
  cell_verts = Geometry.get_faces(topo, D, 0)

  # Fast path: every cell already lists its vertices in increasing global-id
  # order (e.g. simplexified Cartesian models) — π_K = id everywhere.
  all(issorted(cell_verts[c]) for c in 1:length(cell_verts)) && return nothing

  # One reference basis shared by all cells (single-reffe meshes).
  basis = get_prebasis(testitem(cell_reffe))
  rc    = RotationCache(basis)
  cell_ids = IdentityVector(Int32(num_cells(model)))
  cell_change      = lazy_map(_PΛRotationCoBMap(cell_verts, rc, false), cell_ids)
  cell_change_invt = lazy_map(_PΛRotationCoBMap(cell_verts, rc, true),  cell_ids)
  (cell_change, cell_change_invt)
end
