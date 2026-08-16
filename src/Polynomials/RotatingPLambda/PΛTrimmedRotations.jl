# RotatingPLambda/PΛTrimmedRotations.jl
#
# Trimmed (P_r⁻Λ¹) closed form for the generic rotation API of
# RotatingPLambda/PΛRotations.jl: the per-entry `_rotate_basis_function`
# method for `TrimmedPΛBasis` entries `(F, (e1,e2), α)`, plus the pair
# sign/sort helpers (triangle identity, shift operators ρ, anchored pairs).

"""
    trimmed_pair_sign(e1, e2) -> Float64

ε(e1,e2) := 2·[e1<e2] − 1.
"""
trimmed_pair_sign(e1::Int, e2::Int) = e1 < e2 ? 1.0 : -1.0

"""
    trimmed_pair_sort(e1, e2) -> (Int,Int)

e↑ := (min(e1,e2), max(e1,e2)).
"""
trimmed_pair_sort(e1::Int, e2::Int) = minmax(e1, e2)

# Trimmed closed form: single signed term (ε) when the rotated sorted pair is
# anchored (min(π(e)) = min(π(F))); otherwise — uniformly in α, since ϕ does
# not depend on α — the two-term shift expansion
#   ε · ( w(λ;π(F),(m,eπ2),ρ(m,eπ1)∘π(α)) − w(λ;π(F),(m,eπ1),ρ(m,eπ2)∘π(α)) ).
function _rotate_basis_function(
  entries::Vector{Tuple{Vector{Int},Tuple{Int,Int},Vector{Int}}}, idx, w::Int, π::Vector{Int})

  F, e, α  = entries[w]
  e1, e2   = e
  πF       = rotate_face_set(F, π)
  πe1,πe2  = π[e1], π[e2]
  πα       = rotate_multiindex(α, π)
  ε        = trimmed_pair_sign(πe1, πe2)
  eπ1,eπ2  = trimmed_pair_sort(πe1, πe2)
  m      = minimum(πF)

  if m != eπ1   # hit: min(π(e)) ≠ min(π(F))
    # m < eπ1 ⟹ m ∉ {eπ1,eπ2} ⟹ m ∈ supp(π(α)) by the covering condition
    @assert πα[m] > 0 "invalid trimmed entry: hit at α with α[min π(F)] = 0"
    α1 = copy(πα); α1[m] -= 1; α1[eπ1] += 1   # ρ(m,eπ1)∘π(α)
    α2 = copy(πα); α2[m] -= 1; α2[eπ2] += 1   # ρ(m,eπ2)∘π(α)
    [(ε, idx[(πF, (m,eπ2), α1)]), (-ε, idx[(πF, (m,eπ1), α2)])]
  else
    [(ε, idx[(πF, (eπ1,eπ2), πα)])]
  end
end
