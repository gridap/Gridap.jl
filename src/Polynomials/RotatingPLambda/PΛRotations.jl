# RotatingPLambda/PΛRotations.jl
#
# Change of basis induced on the rotating bases by a vertex relabeling
# (rotation) π : ξ → λ, λ = π(ξ). This file implements the closed-form index
# calculus.
#
# No new basis type is introduced: rotating a basis just expresses each of its
# basis functions, pulled back under π⁻¹, as a (signed) linear combination of
# the basis functions of the *same* basis instance, now read with the
# relabeled vertices λ = π(ξ).
#
# The generic API (`bubble_entries`, `bubble_index`, `rotate_basis_function`,
# `rotation_change_of_basis`) is shared by the untrimmed `BarycentricPΛBasis`
# restricted to 1-forms (this file) and the trimmed `BarycentricPmΛBasis`
# (RotatingPLambda/PΛTrimmedRotations.jl); only the per-entry closed form
# `_rotate_basis_function` differs, dispatched on the entry tuple type
# (`k::Int` untrimmed, `e::Tuple{Int,Int}` trimmed).

const _PΛBases = Union{BarycentricPΛBasis,BarycentricPmΛBasis}

##########################
# Vertex/multi-index maps #
##########################

"""
    rotate_face_set(F::Vector{Int}, π::Vector{Int}) -> Vector{Int}

π(F) as a sorted vertex-set key, `sort(π[F])`. The rotating bases never track
face orientation — their bubbles are built from `combinations`, which already
returns canonical sorted vertex sets — so sorting the rotated image loses
nothing this module relies on. Do not reuse this for orientation-sensitive
code, where a face's vertex order encodes its orientation.
"""
rotate_face_set(F::Vector{Int}, π::Vector{Int}) = sort(π[F])

"""
    rotate_multiindex(α, π) -> Vector{Int}

π(α), per `π(α)_i = α_{π⁻¹(i)}`. The only inverse used anywhere in this file:
the exact combinatorial `invperm`, never a numerical/matrix inverse. To apply
π⁻¹ instead, call this (or anything downstream of it) with `invperm(π)`.
"""
rotate_multiindex(α::AbstractVector{Int}, π::Vector{Int}) = α[invperm(π)]

#################################
# Bubble (F,k,α) ↔ index lookup #
#################################

_bubble_entry_type(::BarycentricPΛBasis)  = Tuple{Vector{Int},Int,Vector{Int}}
_bubble_entry_type(::BarycentricPmΛBasis) = Tuple{Vector{Int},Tuple{Int,Int},Vector{Int}}

# A Barycentric bubble function stores its direction form as the index set J.
# For 1-forms that is the single vertex the untrimmed closed form calls k, resp.
# the pair (e1,e2) spanning the Whitney form of the trimmed one.
_bubble_key(::BarycentricPΛBasis,  J) = J[1]
_bubble_key(::BarycentricPmΛBasis, J) = (J[1], J[2])

# The untrimmed law holds for either flavor. The trimmed one expands a hit into
# two terms whose multi-indices differ from α, so the multinomial coefficient of
# Bα does not cancel between them: it holds for the bare monomials λ^α only.
_rotation_flavors(::BarycentricPΛBasis)  = (:AFW, :BMM)
_rotation_flavors(::BarycentricPmΛBasis) = (:BMM,)

"""
    bubble_entries(b) -> Vector{Tuple{Vector{Int},kT,Vector{Int}}}

`(F,k,α)` for each basis function `w`, indexed by `w`. `k` is a single vertex
for the untrimmed basis and a pair `(e1,e2)` for the trimmed one.
"""
function bubble_entries(b::_BaryPΛBasis)
  @check isone(b.k) "The rotation API is only defined for 1-forms, got k=$(b.k)"
  @check b.flavor in _rotation_flavors(b) """
    The rotation API of $(nameof(typeof(b))) is only defined for flavor in \
    $(_rotation_flavors(b)), got $(b.flavor)"""
  entries = Vector{_bubble_entry_type(b)}(undef, length(b))
  for (F, bubble_functions) in get_bubbles(b), (w, α, _, J) in bubble_functions
    entries[w] = (F, _bubble_key(b, J), α)
  end
  entries
end

"""
    bubble_index(b) -> Dict{entry,Int}

Inverse of [`bubble_entries`](@ref): `(F,k,α) → w`.
"""
bubble_index(b::_PΛBases) = Dict(e => w for (w, e) in enumerate(bubble_entries(b)))

#####################
# Basis-level rotate #
#####################

"""
    rotate_basis_function(b, w::Int, π::Vector{Int}) -> Vector{Tuple{Float64,Int}}

Pullback under π⁻¹ of the `w`-th basis function of `b`, expressed as a list of
`(coefficient, target index)` pairs in the same basis `b`, now read with the
relabeled vertices λ = π(ξ). Implements the two cases of the pullback theorems
(single term with coefficient +1, resp. ±ε trimmed; and the filter-hit
expansion with coefficients −1, resp. ±ε trimmed).

For repeated calls (e.g. assembling [`rotation_change_of_basis`](@ref)),
precompute `entries = bubble_entries(b)` and `idx = bubble_index(b)` once and
call the internal `Gridap.Polynomials._rotate_basis_function(entries, idx, w, π)`.
"""
function rotate_basis_function(b::_PΛBases, w::Int, π::Vector{Int})
  entries = bubble_entries(b)
  idx     = bubble_index(b)
  _rotate_basis_function(entries, idx, w, π)
end

# Untrimmed closed form: single term (+1) unless supp(α) = F and
# π(k) = min(π(F)), in which case the resummation identity gives |F|−1 terms
# with coefficient −1.
function _rotate_basis_function(entries::Vector{Tuple{Vector{Int},Int,Vector{Int}}}, idx, w::Int, π::Vector{Int})
  F, k, α = entries[w]
  πF = rotate_face_set(F, π)
  πk = π[k]
  πα = rotate_multiindex(α, π)
  full_support = α[k] > 0   # [s(α)] == F  ⟺  k ∈ [s(α)]  ⟺  α[k] > 0
  if full_support && πk == minimum(πF)
    [(-1.0, idx[(πF, vi, πα)]) for vi in πF if vi != πk]
  else
    [(1.0, idx[(πF, πk, πα)])]
  end
end

#####################
# Cached rotations   #
#####################

"""
    RotationCache(b)

Precomputes the bubble entry table and its index once, and memoises, per
permutation `π`, the signed index map of the rotation:
[`rotation_map`](@ref)`(rc, π)[w]` is the `Vector{Tuple{Float64,Int}}` of
`(coefficient, target index)` pairs of `rotate_basis_function(b, w, π)`.

At most `(D+1)!` distinct permutations exist, each map is computed once and
shared thereafter, and the transformation matrix is never assembled (only the
one-or-few-entry rows are stored). This is the object mesh-level code should
hold, keyed by the per-cell permutation.
"""
struct RotationCache{B<:_PΛBases,E}
  basis   :: B
  entries :: Vector{E}
  idx     :: Dict{E,Int}
  maps    :: Dict{Vector{Int},Vector{Vector{Tuple{Float64,Int}}}}
end

function RotationCache(b::_PΛBases)
  entries = bubble_entries(b)
  idx     = Dict(e => w for (w, e) in enumerate(entries))
  RotationCache(b, entries, idx,
                Dict{Vector{Int},Vector{Vector{Tuple{Float64,Int}}}}())
end

"""
    rotation_map(rc::RotationCache, π::Vector{Int}) -> Vector{Vector{Tuple{Float64,Int}}}

The signed index map of the rotation `π`, one row per basis function, computed
on first request and memoised in `rc`. Rows alias the cache; do not mutate.
"""
function rotation_map(rc::RotationCache, π::Vector{Int})
  m = get(rc.maps, π, nothing)
  m === nothing || return m
  m = [_rotate_basis_function(rc.entries, rc.idx, w, π) for w in 1:length(rc.entries)]
  rc.maps[copy(π)] = m
  m
end

"""
    rotation_change_of_basis(b, π::Vector{Int}) -> Matrix{Float64}

Dense change-of-basis matrix `C` such that the `w`-th basis function of `b`,
pulled back under π⁻¹, equals `∑_w′ C[w,w′] · (w′-th basis function of b, read
with λ = π(ξ))`. Built entirely from [`rotate_basis_function`](@ref); never
calls `inv` on a matrix. To pull back by π⁻¹ instead, call this with
`invperm(π)`.
"""
function rotation_change_of_basis(b::_PΛBases, π::Vector{Int})
  n = length(b)
  entries = bubble_entries(b)
  idx     = bubble_index(b)
  C = zeros(Float64, n, n)
  for w in 1:n
    for (c, w′) in _rotate_basis_function(entries, idx, w, π)
      C[w, w′] += c
    end
  end
  C
end

# ─────────────────────────────────────────────────────────────────────────────
# Change of basis of the "virtually sorted" cell at a vertex permutation π
# ─────────────────────────────────────────────────────────────────────────────
#
# CONVENTION (pinned empirically over all relative vertex orderings of two-
# triangle and two-tet meshes — the competing candidates C(π)ᵀ, C(invperm(π))
# and every un-relabelled variant fail the tangential jump / interpolation
# round-trip tests at O(1)):
#
# π maps sorted position → local index.  The virtually sorted geometric map is
# F̃ = F ∘ A⁻¹ with A = A_{invperm(π)} (reference automorphism
# V_j ↦ V_{invperm(π)[j]}), hence the conforming basis is
#
#   φ̃_μ = (F̃⁻¹)^* ŵ_μ = Σ_ν C(π)[μ,ν] · (F⁻¹)^* ŵ_ν,
#
# with C = rotation_change_of_basis(b, π).
#
# RELABELLING: φ̃_μ is supported on the LOCAL face π(F(μ)), while the conformity
# machinery glues local index μ by the local face F(μ).  The shape function
# stored at μ must therefore be φ̃_{s(μ)}, with s the face-wise order-preserving
# reindexing sending the i-th bubble of face F to the i-th bubble of
# sort(invperm(π)(F)).  Every cell adjacent to a face then enumerates that
# face's functions in the bubble order of its sorted global vertex ids — which
# is what makes the identity face-own-dof permutations correct for every pindex.
#
# In the linear_combination convention (out[j] = Σ_i values[i,j]·in[i]):
#
#   M[i,j]    = C(π)[s(j), i]           (basis side)
#   Minv[i,j] = C(invperm(π))[i, s(j)]  (dof side, restores duality)
#
# using C(π)·C(invperm(π)) = I.  Both are read off the memoised sparse
# rotation_map rows — no matrix inversion.
function compute_pλ_change(rc, π)
  invπ = invperm(π)
  n    = length(rc.entries)

  # s(μ): face-wise order-preserving reindexing local entry → virtual entry,
  # mapping the i-th bubble of face F to the i-th bubble of sort(invperm(π)(F)).
  face_to_ws = Dict(F => [bf[1] for bf in bfs] for (F, bfs) in get_bubbles(rc.basis))
  s = Vector{Int}(undef, n)
  for (F, ws) in face_to_ws
    wsG = face_to_ws[sort(invπ[F])]
    for i in eachindex(ws)
      s[ws[i]] = wsG[i]
    end
  end
  sinv = invperm(s)

  M = zeros(Float64, n, n)
  rows = rotation_map(rc, π)
  for μ in 1:n, (c, w′) in rows[s[μ]]
    M[w′, μ] += c
  end

  Minv = zeros(Float64, n, n)
  rows_inv = rotation_map(rc, invπ)
  for w in 1:n, (c, w′) in rows_inv[w]
    Minv[w, sinv[w′]] += c
  end

  return M, Minv
end
