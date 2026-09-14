# Tensor operations involving proper tensor calculus and differential geometry operations.
#
# This includes pointwise algebra such as
# (∧, ι, ⋆, ♭/♯, koszul, …).

# ============================================================
# Generate-time expression builders
# ============================================================

# `sgn * expr`, with the sign folded in rather than emitted as a factor.
_signed(sgn::Int, expr) = sgn > 0 ? expr : :(-$expr)

# Sum of `terms`. `Expr(:call, :+)` with no term is `+()`, which errors, so an
# empty sum is emitted as an explicit `zero(T)`.
_sum_expr(terms, T) = isempty(terms) ? :(zero($T)) : Expr(:call, :+, terms...)

# Leibniz expansion of the determinant of the `K`×`K` matrix whose `(a,b)` entry
# is the expression `entry(a,b)`. `K ≥ 1`.
function _det_expr(entry, K::Int)
  terms = (_signed(levicivita(p), Expr(:call, :*, (entry(a, p[a]) for a in 1:K)...))
           for p in permutations(1:K))
  Expr(:call, :+, terms...)
end

# Determinant of the `K`×`K` submatrix of the second order tensor named by `A`,
# with rows `I` and columns `J`.
#
# The expansion has `K!` terms per call site. Past `K = 4` that makes compilation
# the dominant cost, so the emitted code calls `_minor_k` and loops instead.
function _minor_expr(A::Symbol, I, J)
  K = length(I)
  K == 0 && return :(one(eltype($A)))
  K > 4  && return :(_minor_k($A, $(Tuple(I)), $(Tuple(J))))
  _det_expr((a, b) -> :($A[$(I[a]), $(J[b])]), K)
end

# Determinant of the K×K submatrix of `A` with rows `I` and columns `J`.
function _minor_k(A, I::NTuple{K,Int}, J::NTuple{K,Int}) where K
  K == 0 && return one(eltype(A))
  det(SMatrix{K,K}(ntuple(n -> A[I[(n-1)%K+1], J[(n-1)÷K+1]], Val(K*K))))
end

# ============================================================
# Alternation map 𝔸
# ============================================================

"""
    alternation(a::MultiValue{NTuple{K,D}}) -> SkewSymKTensorValue{V,K,D}
    𝔸(a)

Alternating part `𝔸(a)` of the order `K` tensor `a` whose `K` indices share
the dimension `D` of variance `V=typeof(D)`, defined by

``\bb{A}(a)_I = \\frac{1}{K!} Σ_σ sign(σ) a_{I∘σ}``.

It is a skew-symmetric tensor. `𝔸` is the identity on such tensors,
so an operation defined through `𝔸` agrees with its `SkewSymKTensorValue`
method on the tensors both accept.
"""
@generated function alternation(a::MultiValue{NTuple{K,Da},T,K}) where {K,Da,T}
  D = Int(Da)
  L = binomial(D, K)
  V = typeof(Da)
  iszero(L) && return :( zero(SkewSymKTensorValue{$V,$K,$D,$T,0}) )

  f = factorial(K)
  comps = Expr[]
  for I in sorted_combinations(D, K)
    terms = [_signed(levicivita(p), :( a[$(I[p]...)] )) for p in permutations(1:K)]
    s = _sum_expr(terms, T)
    push!(comps, isone(f) ? s : :( $s / $f ))
  end

  quote
    SkewSymKTensorValue{$V,$K,$D}($(comps...))
  end
end

alternation(a::SkewSymKTensorValue) = a
const 𝔸 = alternation

# ============================================================
# Exterior (wedge) product ∧
# ============================================================

"""
    ∧(a::ExteriorFormValue{K1,D},   b::ExteriorFormValue{K2,D})
    ∧(a::ExteriorVectorValue{K1,D}, b::ExteriorVectorValue{K2,D})

Pointwise exterior (wedge) product, a `K1+K2`-form or -vector.
"""
@generated function ∧(a::SkewSymKTensorValue{V,K1,D,T1},
                      b::SkewSymKTensorValue{V,K2,D,T2}) where {V,K1,K2,D,T1,T2}
  D = Int(D)
  K = K1 + K2
  L = binomial(D, K)
  T = Base.promote_op(*, T1, T2)

  iszero(L) && return :( zero(SkewSymKTensorValue{V,$K,D,$T}) )

  # (a ∧ b)_I = Σ ε(I₁,I₂) a_{I₁} b_{I₂}, summed over the ways of splitting the
  # sorted I into a K1-subset I₁ and a K2-subset I₂. Overlapping splits have
  # ε = 0 and contribute nothing.
  terms = [Expr[] for _ in 1:L]
  for (n1, I1) in enumerate(sorted_combinations(D, K1)),
      (n2, I2) in enumerate(sorted_combinations(D, K2))
    s = sorting_sign(I1..., I2...)
    s == 0 && continue
    l = combination_index(sort([I1..., I2...]), D)
    push!(terms[l], _signed(s, :(indep_comp_getindex(a,$n1) * indep_comp_getindex(b,$n2))))
  end
  comps = [_sum_expr(t, T) for t in terms]

  quote
    SkewSymKTensorValue{V,$K,$D}($(comps...))
  end
end

# ============================================================
# Interior product ι_v(ω): contracts ω(v, -)
# ι_v : Ωᴷ(D) → Ωᴷ⁻¹(D),   (ι_v ω)(w₂,…,wₖ) = ω(v, w₂,…,wₖ)
# ============================================================

"""
    interior_product(v::VectorValue, ω::ExteriorFormValue)
    interior_product(v, ω::SkewSymKTensorValue)
    ι(v, ω)

Interior product (contraction) `ι_v ω`, of a an alternating `K`-tensor `ω` and
(co)vector `v`. It is the `K-1`-alternating tensor:
    (ι_v ω)(w₂,…,wₖ) = ω(v, w₂,…,wₖ).

If `v` is not a `VectorValue`, the index of `v` should be contractable with
that of `ω`, else, `v` is assumed contravariant.
"""
@generated function interior_product(v::MultiValue{Tuple{Dv},Tv}, ω::SkewSymKTensorValue{V,K,Dw,Tw}) where {V,K,Dv,Dw,Tv,Tw}
  K >= 1 || return :(@unreachable "interior product requires K ≥ 1")
  Dv == Dw || return :(throw(DimensionMismatch()))
  can_contract(Dv, Dw) || return :(throw(ArgumentError("interior_product undefined for size $(size(v)) vector and alternating tensor with variance $V, indices cannot contract")))

  D = Int(Dv)
  Km1 = K - 1
  T = Base.promote_op(*, Tv, Tw)
  iszero(binomial(D, K)) && return :( zero(SkewSymKTensorValue{V,$Km1,$Dw,$T}) )

  L = binomial(D, Km1)

  # (ι_v ω)_{I∖I[j]} = Σ_j (-1)^{j-1} v^{I[j]} ω_I. Dropping one entry of the
  # sorted I leaves it sorted, so the output slot is its combination index.
  terms = [Expr[] for _ in 1:L]
  for (n, I) in enumerate(sorted_combinations(D, K)), j in 1:K
    l = combination_index(deleteat!(copy(I), j), D)
    push!(terms[l], _signed(iseven(j-1) ? 1 : -1,
                            :(v[$(I[j])] * indep_comp_getindex(ω,$n))))
  end
  comps = [_sum_expr(t, T) for t in terms]

  :( SkewSymKTensorValue{V,$Km1,$Dw}(($(comps...),)) )
end

# VectorValue is assumed contravariant for the interior product
function interior_product(v::VectorValue{D}, ω::ExteriorFormValue) where D
  contra_v = ContraVectorValue{D}(Tuple(v))
  interior_product(contra_v, ω)
end
const ι = interior_product

# ============================================================
# Hodge star ⋆ (flat Euclidean metric, g = Iᴅ, so √|det(g)| = 1)
#
# (⋆ω)_J = ∑_I  ε_{IJ}  ω_I
# where I ranges over K-combinations, J over (D-K)-combinations,
# and ε_{IJ} = levicivita of the permutation (I,J) → sorted (1…D).
# Only the complementary pair J = Iᶜ contributes, since ε_{IJ} vanishes as soon
# as I and J share an index.
# ============================================================

"""
    _hodge_star_signs(K, D) -> Vector{Int}

Levi-Civita signs of the flat Hodge star, entry `[n] = ε_{I Iᶜ}` for `I` the
n-th `K`-combination of `1:D`, so that with `L = binomial(D,K)`

    (⋆ω)_{L+1-n} = _hodge_star_signs(K,D)[n] * ω_n.

Only the complementary pair contributes, and its output slot is `L+1-n` because
complementation reverses the lexicographic order of the combinations. The sign
is `(-1)^(ΣI - K(K+1)/2)`, since the sorted `I` followed by its sorted
complement has `Σₐ (I[a] - a)` inversions.
"""
function _hodge_star_signs(K::Int, D::Int)

  # Helper: advance `I` to the next K-combination of 1:D in lexicographic
  # order, leaving the last combination unchanged.
  function _next_combination!(I::Vector{Int}, D::Int)
    K = length(I)
    j = K
    while j >= 1 && I[j] == D - K + j
      j -= 1
    end
    j == 0 && return I
    I[j] += 1
    for l in j+1:K
      I[l] = I[l-1] + 1
    end
    I
  end

  L = binomial(D, K)
  signs = Vector{Int}(undef, L)
  I = collect(1:K)              # first K-combination in lexicographic order
  offset = (K*(K+1)) ÷ 2
  for n in 1:L
    signs[n] = iseven(sum(I) - offset) ? 1 : -1
    _next_combination!(I, D)
  end
  signs
end


"""
    flat_hodge_star(ω::ExteriorFormValue{K,D})
    flat_hodge_star(ω::MultiValue{NTuple{K,D}})
    ⋆ᴱ(ω)

Hodge star `⋆ω` in Euclidean space, that is [`hodge_star`](@ref) with the
identity metric, an alternating `(D-K)`-tensor with same variance as `ω`.
"""
@generated function flat_hodge_star(ω::SkewSymKTensorValue{V,K,Dw}) where {V,K,Dw}
  D = Int(Dw)
  K <= D || return :(@unreachable $("hodge star requires K ≤ D, got K = $K and D = $D"))

  Kc = D - K
  L  = binomial(D, K)   # = binomial(D, Kc), the two are matched by complementation
  s  = _hodge_star_signs(K, D)
  comps = [_signed(s[L+1-m], :(indep_comp_getindex(ω,$(L+1-m)))) for m in 1:L]
  :( SkewSymKTensorValue{V,$Kc,$Dw}($(comps...)) )
end
flat_hodge_star(a::MultiValue{NTuple{K,Dw}}) where {K,Dw} = flat_hodge_star(𝔸(a))
const ⋆ᴱ = flat_hodge_star

# ============================================================
# Generic functions extended elsewhere
#
# Methods for Symbolics.Num coefficients are provided by the GridapSymbolicsExt
# package extension.
# ============================================================

"""
    symbolic_coordinates(D::Integer)

The first `D` global Cartesian symbolic coordinate variables `x¹,…,x⁹`
(at most 9), as a tuple of `Symbolics.Num`. These are the variables
`exterior_derivative` differentiates with respect to for symbolic forms.

Provided by the GridapSymbolicsExt package extension when Symbolics is loaded.
"""
function symbolic_coordinates end

# ============================================================
# Riemannian Hodge star  ⋆_g : Ωᴷ → Ωᴰ⁻ᴷ
#
# Algorithm (coordinate formula):
#   Step 1 — raise indices:
#     ω^I = Σ_{J sorted K-index} det([g_inv[I[a],J[b]]]_{a,b}) ω_J
#   Step 2 — Levi-Civita contraction:
#     (⋆_g ω)_L = √|det(g)| Σ_{I sorted} ε(I,L) ω^I
#
# Passing g_inv (the pointwise inverse metric tensor) and sqrt_abs_g
# (the pointwise square root of det g) keeps this purely algebraic.
# ============================================================

"""
    hodge_star(ω::MultiValue{NTuple{K,D}},  g_or_ginv[, sqrt_abs_g])
    ⋆(ω, g_or_ginv[, sqrt_abs_g])

Hodge star `⋆ω` of a `K`-tensor `ω` with variance. It is an alternating
`(D-K)`-tensor of the same variance, a form for a covariant `ω` and a
`(D-K)`-vector for a contravariant one.

`g_or_ginv` is either the metric tensor `g` (covariant) or its inverse `g_inv`
(contravariant); whichever is missing is computed on demand, see
[`parse_metric_data`](@ref). Passing `sqrt_abs_g` = `√|det(g)|` spares computing
a determinant.

A not-alternating `K`-tensor with homogeneous variance is accepted, and
starred through its alternating part.

See also [`flat_hodge_star(ω)`](@ref) for the Euclidean version.
"""
@generated function hodge_star(ω::SkewSymKTensorValue{V,K,Dw,Tw}, g_or_ginv::SymTensorValue{Dg},
                               sqrt_abs_g=nothing) where {V,K,Dg,Dw,Tw}
  V === Int && return :(@unreachable $("hodge_star requires a K-tensor with variance, got $V index"))
  Dg == Dw || return :(throw(DimensionMismatch()))
  D = Int(Dg)
  K <= D || return :(@unreachable $("hodge star requires K ≤ D, got K = $K and D = $D"))

  Kc  = D - K
  L   = binomial(D, K)   # = binomial(D, Kc), the two are matched by complementation
  c_K = sorted_combinations(D, K)

  # Raising the indices of a form pairs with the volume factor √|det(g)|;
  # lowering those of a K-vector pairs with its reciprocal.
  _g_or_ginv = V === Cov ? :g_inv : :g
  _volume    = V === Cov ? :sqrt_abs_g : :(inv(sqrt_abs_g))

  # Step 1: swap (raise or lower) K indices of ω. `K ≤ D` makes the sum non-empty.
  swapped = [Expr(:call, :+, (:( $(_minor_expr(_g_or_ginv, I, J)) * indep_comp_getindex(ω,$m) )
                             for (m, J) in enumerate(c_K))...)
            for I in c_K]

  # Step 2: Levi-Civita contraction, which pairs each I with its complement only
  s = _hodge_star_signs(K, D)
  comps = [_signed(s[L+1-m], :(vol * ω_swapped[$(L+1-m)])) for m in 1:L]

  # The metric data is read into locals once: `MetricData` is immutable, so a
  # field it does not store is recomputed at every access.
  quote
    (; $_g_or_ginv, sqrt_abs_g) = parse_metric_data(g_or_ginv, sqrt_abs_g)
    vol = $_volume
    ω_swapped = ($(swapped...),)
    SkewSymKTensorValue{$V,$Kc,$Dw}(($(comps...),))
  end
end

function hodge_star(a::MultiValue{NTuple{K,D}}, g_or_ginv, sqrt_abs_g=nothing) where {K,D}
  hodge_star(𝔸(a), g_or_ginv, sqrt_abs_g)
end

const ⋆ = hodge_star

# ============================================================
# Musical isomorphisms  ♭ (flat) and ♯ (sharp)
#
# flat  (♭):  ContraVectorValue{D} → CoVectorValue{D}
#   (v♭)_i = Σ_j g_{ij} v^j           (lowers an index)
#
# sharp (♯):  CoVectorValue{D} → ContraVectorValue{D}
#   (ω♯)^i = Σ_j g^{ij} ω_j           (raises an index)
# ============================================================

"""
    flat(v::VectorValue,          g_or_ginv::SymTensorValue)
    flat(v::MultiValue{Tuple{D}}, g_or_ginv::SymTensorValue)
    ♭(v, g_or_ginv)

Musical isomorphism ♭: lower the index of the contravariant vector `v` with the
metric `g`, resulting into a covector (1-form).

`g_or_ginv` is either the metric tensor `g` (covariant) or its inverse `g_inv`
(contravariant), see [`parse_metric_data`](@ref).
"""
@generated function flat(v::MultiValue{Tuple{Dv},Tv}, g_or_ginv::SymTensorValue{Dg,Tg}) where {Dv,Dg,Tv,Tg}
  Dg == Dv || return :(throw(DimensionMismatch()))
  Dv isa Con || return :(throw(ArgumentError("flat undefined for 1-tensor that is not contravariant, got $Dv")))
  D = Int(Dv)
  Dcov = inverse_variance(Dv)
  iszero(D) && return :( zero(CoVectorValue{$Dcov,$(Base.promote_op(*,Tg,Tv))}) )

  comps = [Expr(:call, :+, (:(g[$i,$j] * v[$j]) for j in 1:D)...) for i in 1:D]
  quote
    (; g) = parse_metric_data(g_or_ginv)
    CoVectorValue{$Dcov}(($(comps...),))
  end
end
function flat(v::VectorValue{D}, g_or_ginv::SymTensorValue) where D
  contra_v = ContraVectorValue{D}(Tuple(v))
  flat(contra_v, g_or_ginv)
end

const ♭ = flat

"""
    sharp(ω, g_or_ginv::SymTensorValue{D})
    ♯(ω, g_or_ginv)

Musical isomorphism `♯`: raise the index of the covector (1-form) `ω` with the
inverse metric `g_inv`, resulting into a contravariant vector.

`g_or_ginv` is either the metric tensor `g` (covariant) or its inverse `g_inv`
(contravariant), see [`parse_metric_data`](@ref).
"""
@generated function sharp(ω::MultiValue{Tuple{Dw},Tw}, g_or_ginv::SymTensorValue{Dg,Tg}) where {Dg,Dw,Tw,Tg}
  Dg == Dw || return :(throw(DimensionMismatch()))
  Dw isa Cov || return :(throw(ArgumentError("sharp undefined for 1-tensor that is not covariant, got $Dw")))
  D = Int(Dw)
  Dcon = inverse_variance(Dw)
  iszero(D) && return :( zero(ContraVectorValue{$Dcon,$(Base.promote_op(*,Tg,Tw))}) )

  comps = [Expr(:call, :+, (:(g_inv[$i,$j] * indep_comp_getindex(ω,$j)) for j in 1:D)...) for i in 1:D]
  quote
    (; g_inv) = parse_metric_data(g_or_ginv)
    ContraVectorValue{$Dcon}(($(comps...),))
  end
end

const ♯ = sharp

"""
    apply_tensor(ω::ExteriorFormValue{K,D,T}, v₁, …, vₖ) → scalar

Evaluate the K-form ω on K tangent contravariant vectors v₁,…,vₖ ∈ ℝᴰ:
  ω(v₁,…,vₖ) = Σ_{|I|=K} ω_I · det([vₐ[I[b]]]_{a,b=1}^K)

Each vector must be a size `(D, )` tensor. Behaves like `getindex` for scalar `K=0`.
"""
@generated function apply_tensor(ω::SkewSymKTensorValue{V,K,Dw,T},
                                 vs::Vararg{MultiValue{Tuple{D}},K}) where {V,K,D,Dw,T}
  K == 0 && return :(ω[1])
  D == Dw || return :(throw(DimensionMismatch()))
  can_contract(D,Dw) || return :(throw(ArgumentError("tensor $(ω) cannot be applied to  $(vs), indices cannot contract")))

  # trivial Λᴷ for K > D
  iszero(binomial(D, K)) &&
    return Expr(:call, :*, :(zero($T)), (:(zero(eltype(vs[$a]))) for a in 1:K)...)

  terms = [:( indep_comp_getindex(ω,$n) * $(_det_expr((a,b) -> :(vs[$a][$(I[b])]), K)) )
    for (n, I) in enumerate(sorted_combinations(Int(D), K))]
  Expr(:call, :+, terms...)
end
function apply_tensor(ω::ExteriorFormValue{K}, vs::Vararg{VectorValue{D},K}) where {K,D}
  K == 0 && return ω[1]
  contra_vs = map(vi -> ContraVectorValue{D}(Tuple(vi)), vs)
  apply_tensor(ω, contra_vs...)
end

# TODO: generalize apply_tensor to non skew-symmetric tensor ω

apply_tensor(ω::SkewSymKTensorValue{V,K,D}, vs...) where {V,K,D} =
  @unreachable "a K-tensor takes K arguments (subtyping `MultiValue`) of size (D,), got $(typeof.(vs)) for K = $K and D = $D"

# ============================================================
# Pullback  φ*ω : Ωᴷ(Dn) → Ωᴷ(Dm)
#
# Given a map φ : ℝᴰᵐ → ℝᴰⁿ from a chart to an ambient space, its transposed
# Jacobian Jt = ∇φ, Jt[j,i] = ∂φⁱ/∂ξʲ, and a K-form ω on ℝᴰⁿ:
#
#   (φ*ω)_J(ξ) = Σ_I  ω_I(φ(ξ)) · det(Jt[J, I])
#
# where I is a sorted K-multi-index into {1:Dn} (ambient) and
# J is a sorted K-multi-index into {1:Dm} (chart).
# ============================================================

"""
    pullback(ω::ExteriorFormValue{K,Dn}, Jt::TensorValue{Dm,Dn})

Pointwise pullback `φ^*ω` of a K-form under a map `φ : ℝᴰᵐ → ℝᴰⁿ` with
transposed Jacobian `Jt = ∇φ`, that is `Jt[j,i] = ∂φⁱ/∂ξʲ`.

The chart index `j` of `Jt` is covariant and its ambient index `i`
contravariant, so `Jt` has size `(Cov(Dm), Con(Dn))`: the ambient index is the
one contracted with `ω`, and the result is an `ExteriorFormValue{K,Cov(Dm)}`.

A Field-level method for lazy pullbacks is provided in `Gridap.Fields`.
"""
@generated function pullback(ω::ExteriorFormValue{K,Dw,Tw}, Jt::TensorValue{Dm,Dn,Tj}) where {K,Dw,Dm,Dn,Tw,Tj}
  Dm isa Cov || return :(throw(ArgumentError(
    "the chart index of Jt is covariant, got a transposed jacobian of size $(size(Jt))"
  )))
  can_contract(Dn, Dw) || return :(throw(DimensionMismatch(
    "cannot pullback a form of dimension $Dw by a transposed jacobian of ambient dimension $Dn, indices cannot contract"
  )))
  T = Base.promote_op(*, Tw, Tj)

  m, n = Int(Dm), Int(Dn)
  (iszero(binomial(n, K)) || iszero(binomial(m, K))) &&
    return :( zero(ExteriorFormValue{$K,$Dm,$T}) )

  c_I = sorted_combinations(n, K)   # K-combins in ambient (target) space
  c_J = sorted_combinations(m, K)   # K-combins in chart (origin) space

  # det([∂φ^{I[b]}/∂ξ^{J[a]}]) is the minor of Jt on rows J and columns I
  comps = [_sum_expr([:( $(_minor_expr(:Jt, J, I)) * indep_comp_getindex(ω,$l) )
                      for (l, I) in enumerate(c_I)], T)
           for J in c_J]

  :( ExteriorFormValue{$K,$Dm}(($(comps...),)) )
end
# TODO: generalize the pullback to tensor of type (r,s,0)
# - mathematically express the general pullback using contractions with J and Jt_pinv, add Jt_pinv=pinvJt(Jt) as optional argument
# - implement using tensor_contraction

# TODO implement the pushforward of k-forms, this is just taking the inverse of Jt and feeding to pullback

# TODO implement the pushforward of type (r,s,0) tensors, this is just taking the inverse of Jt and feeding to pullback

# ============================================================
# Koszul operator  κ_x(ω) = ι_x(ω)  (value level)
#
# κ : Λᴷ → Λᴷ⁻¹,  κ_x(ω) = interior_product(x, ω)
#
# Homotopy identity on Pᵣ Λᴷ: dκ + κd = (r+K) id
# This is the fundamental algebraic fact behind trimmed polynomial
# form spaces: P⁻ᵣΛᴷ = Pᵣ₋₁Λᴷ ⊕ κ(Hᵣ₋₁Λᴷ⁺¹)
# ============================================================

"""
    koszul(x::VectorValue, ω::ExteriorFormValue)

Koszul differential ``κ_x(ω) = ι_x(ω)``, the interior product with the position
vector `x`. A Field-level method is provided in `Gridap.Fields`.
"""
function koszul(x::VectorValue{D}, ω::ExteriorFormValue{K,Dw}) where {K,D,Dw}
  interior_product(x, ω)
end

# ============================================================
# Inner product of K-forms
#
# (α|β) = (1/K!) α_{i₁…i_K} g^{i₁j₁} ⋯ g^{i_Kj_K} β_{j₁…j_K}
#       = Σ_{|I|=|J|=K} α_I β_J det([g⁻¹[I[a],J[b]]]_{a,b})
#
# the second sum running over sorted multi-indices, which reduces to
# Σ_I α_I β_I for a flat metric. This is the inner product the Hodge star is
# built on,
#   α ∧ ⋆_g β = (α|β) μ = (α|β) √|det(g)| dx¹ ∧ … ∧ dx^D.
#
# `flat_inner` is the tensor inner product with the indices of β raised by the
# identity metric, such that flat_inner(α,β) = K!(α|β).
# ============================================================

"""
    form_inner(ω::ExteriorFormValue{K,D}, η::ExteriorFormValue{K,D}, g_or_ginv::SymTensorValue{D})
    ⨟(ω,η)
    ω ⨟ η

Inner product ``(ω|η)`` of two K-forms, the scalar such that

``ω ∧ ⋆_g η = √|det(g)| (ω|η) dx¹∧…∧dx^D``.

`g_or_ginv` is either the metric tensor `g` (covariant) or its inverse `g_inv`
(contravariant), see [`parse_metric_data`](@ref). See
[`flat_form_inner`](@ref) for the identity metric.
"""
@generated function form_inner(a::ExteriorFormValue{K,Dw,Ta}, b::ExteriorFormValue{K,Dw,Tb},
                               g_or_ginv::SymTensorValue{Dg,Tg}) where {K,Dg,Dw,Ta,Tb,Tg}
  Dg == Dw || return :(throw(DimensionMismatch()))
  D = Int(Dg)
  T = Base.promote_op(*, Ta, Tb, Tg)
  c_K = sorted_combinations(D, K)
  # raising β_I to β^I as `hodge_star` does
  terms = [:( indep_comp_getindex(a,$m) *
              $(_sum_expr([:( $(_minor_expr(:g_inv, I, J)) * indep_comp_getindex(b,$n) )
                           for (n, J) in enumerate(c_K)], T)) )
           for (m, I) in enumerate(c_K)]
  # (α|β) = Σ_I α_I β^I
  quote
    (; g_inv) = parse_metric_data(g_or_ginv)
    $(_sum_expr(terms, T))
  end
end
# TODO: same as hodge_star, generalize to SkewSymKTensorValue,
const ⨟ = form_inner

"""
    flat_form_inner(ω::ExteriorFormValue{K,D},   η::ExteriorFormValue{K,D})
    flat_form_inner(ω::ExteriorVectorValue{K,D}, η::ExteriorVectorValue{K,D})
    ⨟ᴱ(ω,η)
    ω ⨟ᴱ η

Inner product ``(ω|η)`` of two `K`-forms or -vectors in Euclidean space, that is
[`form_inner`](@ref) with the identity (inverse) metric.

This differs from the tensor inner product by a factor `K!`, see [`flat_inner`](@ref).
"""
@generated function flat_form_inner(
  a::SkewSymKTensorValue{V,K,D,Ta}, b::SkewSymKTensorValue{V,K,D,Tb}) where {V,K,D,Ta,Tb}

  T = Base.promote_op(*,Ta,Tb)
  L = num_indep_components(a)
  _sum_expr([:(indep_comp_getindex(a,$i) * indep_comp_getindex(b,$i)) for i in 1:L], T)
end
const ⨟ᴱ = flat_form_inner

"""
    flat_inner(ω::MultiValue{S}, η::MultiValue{S})
    ⊙ᴱ(ω,η)
    ω ⊙ᴱ η

Tensor inner product in Euclidean space, that is `inner`/`⊙` even if `ω` and
`η` have variance but with iddentity (inverse) metric.
"""
function flat_inner(a::SkewSymKTensorValue{V,K,D,Ta},
                    b::SkewSymKTensorValue{V,K,D,Tb}) where {V,K,D,Ta,Tb}
  K > D && return zero(Base.promote_op(*,Ta,Tb))
  factorial(K) * flat_form_inner(a, b)
end
const ⊙ᴱ = flat_inner

# TODO: generalize flat_inner to tensor of type (r,s,u), you can e.g. invert
# variace of `a` with a `convert`, and call inner, (or remove the variance of
# both argument by `convert` if that's simpler).


# ============================================================
# Isomorphism: VectorValue{D} ↔ CoVectorValue{D}
#
# In a flat Euclidean space, 1-forms and vectors are identified
# via the standard inner product.  These conversions let you
# take a VectorValue FEFunction (e.g., from a Nédélec FESpace)
# and treat it as a genuine differential 1-form for subsequent
# operations (d, ⋆, ∧, …).
# ============================================================

"""
    grad_to_2form(Jt::TensorValue{D,D})

Exterior derivative of a proxied 1-form `ω::VectorValue` from its gradient `Jt
= ∇ω`:

``(dω)_{a<b} = ∂ω_b/∂x^a − ∂ω_a/∂x^b = Jt[a,b] − Jt[b,a].``

Used to compute the exterior derivative of a vector proxied 1-form,, see `d_1form` in `Gridap.CellData`.
"""
@generated function grad_to_2form(Jt::TensorValue{D,D,T}) where {D,T}
  D isa Int || return :( @unreachable $(
    "grad_to_2form requires un-varianced indices, got a TensorValue{$D,$D}") )
  iszero(binomial(D, 2)) && return :( zero(ExteriorFormValue{2,$D,$(Base.promote_op(-,T,T))}) )

  comps = [:( Jt[$a,$b] - Jt[$b,$a] ) for (a, b) in sorted_combinations(D, 2)]
  :( ExteriorFormValue{2,$D}(($(comps...),)) )
end


"""
    to_1form(v::VectorValue{D})

Conversion of `v` to a 1-form, a `CoVectorValue{D}`, assuming Euclidean space
(the metric tensor is the iddentity matrix).
"""
to_1form(v::VectorValue{D,T}) where {D,T} = CoVectorValue{D,T}(Tuple(v))

"""
    from_1form(ω::CoVectorValue)

Conversion of the 1-form value `ω` to a `VectorValue{D}`, assuming Euclidean space
(the metric tensor is the iddentity matrix).
"""
from_1form(ω::CoVectorValue{D,T}) where {D,T} = VectorValue{Int(D),T}(Tuple(ω))

# General K-form ↔ VectorValue{binomial(D,K)} (component-wise isomorphism)

"""
    to_Kform(v::VectorValue{L}, ::Val{K}, ::Val{D})

Reinterpret the `L = binomial(D,K)` components of `v` as the components of a
`ExteriorFormValue{K,D}` in canonical Cartesian basis `dx^I`.
"""
to_Kform(v::VectorValue{L,T}, ::Val{K}, ::Val{D}) where {L,T,K,D} = ExteriorFormValue{K,D,T}(Tuple(v))

"""
    from_Kform(ω::ExteriorFormValue{K,D})

Reinterpret the components of `ω` as a `VectorValue{binomial(D,K)}`, the inverse
of [`to_Kform`](@ref).
"""
from_Kform(ω::ExteriorFormValue{K,D,T}) where {K,D,T} = VectorValue{binomial(Int(D),K),T}(Tuple(ω))

# ============================================================
# Scalar ↔ 0-form and D-form
# ============================================================

"""
    to_0form(u::Number, ::Val{D})

Convert the scalar `u` to a `ExteriorFormValue{0,D}`.
"""
to_0form(u::_Scalar, ::Val{D}) where D = ExteriorFormValue{0,D}((u,))

"""
    to_Dform(u::Number, ::Val{D})

Convert the scalar `u` as the sole component of the top-degree
`ExteriorFormValue{D,D}`, the inverse of [`vol_coeff`](@ref).
"""
to_Dform(u::_Scalar, ::Val{D}) where D = ExteriorFormValue{D,D}((u,))

# ============================================================
# Volume form scalar extraction
#
# vol_coeff : Λᴰ → ℝ
# A top-degree D-form in D dimensions has exactly one component.
# Extracting it gives the scalar integrand for ∫(ω ∧ ⋆η) dΩ.
# ============================================================

"""
    vol_coeff(ω::ExteriorFormValue{D,D})

The single scalar coefficient of a top-degree D-form in D dimensions.
Enables using Gridap's standard measure for integration of differential forms:
    ∫(vol_coeff(ω ∧ ⋆η)) * dΩ
""" # to integrate a k-vector, the hodge_star (with metric) needs to be applied before vol_coeff.
function vol_coeff(ω::ExteriorFormValue{K,D}) where {K,D}
  K == D || @unreachable "vol_coeff requires a top-degree form, got K = $K and D = $(Int(D))"
  indep_comp_getindex(ω, 1)
end

