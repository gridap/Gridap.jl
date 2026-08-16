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
# `rotation_change_of_basis`) is shared by `RotatingPΛBasis` (untrimmed, this
# file) and `TrimmedPΛBasis` (RotatingPLambda/PΛTrimmedRotations.jl); only the
# per-entry closed form `_rotate_basis_function` differs, dispatched on the
# entry tuple type (`k::Int` untrimmed, `e::Tuple{Int,Int}` trimmed).

const _PΛBases = Union{RotatingPΛBasis,TrimmedPΛBasis}

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

_bubble_entry_type(::RotatingPΛBasis) = Tuple{Vector{Int},Int,Vector{Int}}
_bubble_entry_type(::TrimmedPΛBasis)  = Tuple{Vector{Int},Tuple{Int,Int},Vector{Int}}

"""
    bubble_entries(b) -> Vector{Tuple{Vector{Int},kT,Vector{Int}}}

`(F,k,α)` for each basis function `w`, indexed by `w`. `k` is a single vertex
for `RotatingPΛBasis` and a pair `(e1,e2)` for `TrimmedPΛBasis`.
"""
function bubble_entries(b::_PΛBases)
  entries = Vector{_bubble_entry_type(b)}(undef, length(b))
  for (F, bubble_functions) in b.bubbles, (w, k, α, _) in bubble_functions
    entries[w] = (F, k, α)
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
