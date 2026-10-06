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
    div_op = T <: Integer ? div : / # preserve ints eltype
    push!(comps, isone(f) ? s : :( $(div_op)($(s), $(f)) ))
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
    ∧(a::MultiValue{NTuple{K1,D}},  b::MultiValue{NTuple{K1,D}})
    a ∧ b

Pointwise exterior (wedge) product, a `K1+K2`-form or -vector. It is formaly
defined by

    a ∧ b = (K1+K2)! / (K1!*K2!) * 𝔸(a⊗b)
          = 𝔸(a) ∧ 𝔸(b)

but is computed more efficiently if the arguments are alternating, that is are
either K-forms or -vectors.
"""
function ∧(a, b)
  throw(ArgumentError(
    "wedge product undefined for these sizes, all dimension indices \
    should be iddentical (same number and variance), got sizes $(size(a)) and $(size(b)).
  "))
end

∧(a::MultiValue{NTuple{K1,D}}, b::MultiValue{NTuple{K2,D}}) where {K1,K2,D} = 𝔸(a)∧𝔸(b)
# Next method is needed if K1 or K2=0 and Ds are different, since that can
# dispatch to the generic but not in the altinating tensor implementation
function ∧(a::SkewSymKTensorValue{V1,K1}, b::SkewSymKTensorValue{V2,K2}) where {V1,V2,K1,K2}
  throw(ArgumentError("wedge product undefined for tensors of different variances, got sizes $(size(a)) and $(size(b))."))
end

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
that of `ω`, else, `v` is assumed contravariant so `ω` should be covariant.
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
  for (n, I) in enumerate(sorted_combinations(D, K)),
      (j, l) in enumerate(sub_combinations_ids(I, D))
    push!(terms[l], _signed(iseven(j-1) ? 1 : -1,
                            :(v[$(I[j])] * indep_comp_getindex(ω,$n))))
  end
  comps = [_sum_expr(t, T) for t in terms]

  :( SkewSymKTensorValue{V,$Km1,$Dw}(($(comps...),)) )
end

const ι = interior_product


"""
    symbolic_coordinates(D::Integer)

The first `D` global Cartesian symbolic coordinate variables `x¹,…,x⁹`
(at most 9), as a tuple of `Symbolics.Num`. These are the variables
`exterior_derivative` differentiates with respect to for symbolic forms.

Provided by the GridapSymbolicsExt package extension when Symbolics is loaded.
"""
function symbolic_coordinates end

function _check_metric_dim(g_or_ginv::SymTensorValue{Dg}, D::Int) where Dg
  Int(Dg) == D || throw(DimensionMismatch(
    "the metric tensor has dimension $(Int(Dg)), but the tensor it acts on has dimension $D"))
  nothing
end

# ============================================================
# Riemannian Hodge star  ⋆ : Ωᴷ → Ωᴰ⁻ᴷ
#
# Algorithm (coordinate formula):
#   Step 1 — raise indices:
#     ω^I = Σ_{J sorted K-index} det([g_inv[I[a],J[b]]]_{a,b}) ω_J
#   Step 2 — Levi-Civita contraction:
#     (⋆ ω)_L = √|det(g)| Σ_{I sorted} ε(I,L) ω^I
#
# Passing g_inv (the pointwise inverse metric tensor) and sqrt_abs_g
# (the pointwise square root of det g) keeps this purely algebraic.
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
    hodge_star(ω::MultiValue{NTuple{K,D}}, metric_args...)
    ⋆(ω, metric_args...)

Hodge star `⋆ω` of a `K`-tensor `ω` with variance. It is an alternating
`(D-K)`-tensor of the same variance, a form for a covariant `ω` and a
`(D-K)`-vector for a contravariant one.

`metric_args` contains the metric information, see [`parse_metric_data`](@ref).
`hodge_star` uses `sqrt_abs_g` and the metric with opposite variance to `ω`,
that is `g_inv` for a form and `g` for a `K`-vector.

A not-alternating `K`-tensor with homogeneous variance is accepted, and
starred through its alternating part.

See also [`flat_hodge_star(ω)`](@ref) for the Euclidean version, accessible by
passing `LinearAlgebra.I`.
"""
@generated function hodge_star(ω::SkewSymKTensorValue{V,K,Dw,Tw}, metric_args::Vararg{_MetricArg}) where {V,K,Dw,Tw}
  V === Int && return :(@unreachable $("hodge_star requires a K-tensor with variance, got $V index"))
  D = Int(Dw)
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
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_hodge_star(ω)
    (; $_g_or_ginv, sqrt_abs_g) = md
    _check_metric_dim($_g_or_ginv, $D)
    vol = $_volume
    ω_swapped = ($(swapped...),)
    SkewSymKTensorValue{$V,$Kc,$Dw}(($(comps...),))
  end
end

function hodge_star(a::MultiValue{NTuple{K,D}}, metric_args::Vararg{_MetricArg}) where {K,D}
  hodge_star(𝔸(a), metric_args...)
end

const ⋆ = hodge_star

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
# Musical isomorphisms  ♭ (flat) and ♯ (sharp)
#
# flat  (♭):  ContraVectorValue{D} → CoVectorValue{D}
#   (v♭)_i = Σ_j g_{ij} v^j           (lowers an index)
#
# sharp (♯):  CoVectorValue{D} → ContraVectorValue{D}
#   (ω♯)^i = Σ_j g^{ij} ω_j           (raises an index)
# ============================================================

"""
    flat(v::VectorValue,          metric_args...)
    flat(v::MultiValue{Tuple{D}}, metric_args...)
    ♭(v, metric_args...)

Musical isomorphism ♭: lower the index of the contravariant vector `v` with the
metric `g`, resulting into a covector (1-form).

`metric_args` contains the metric information, see [`parse_metric_data`](@ref).
`flat` uses `g`.

See also [`flat_flat(v)`](@ref) for the Euclidean version, accessible by
passing `LinearAlgebra.I`.
"""
@generated function flat(v::MultiValue{Tuple{Dv},Tv}, metric_args::Vararg{_MetricArg}) where {Dv,Tv}
  Dv isa Con || return :(throw(ArgumentError("flat undefined for 1-tensor that is not contravariant, got $Dv")))
  D = Int(Dv)
  Dcov = inverse_variance(Dv)

  # An empty result still promotes with the scalar type of the metric.
  comps = [Expr(:call, :+, (:(g[$i,$j] * v[$j]) for j in 1:D)...) for i in 1:D]
  result = iszero(D) ? :( zero(CoVectorValue{$Dcov,Base.promote_op(*,eltype(g),$Tv)}) ) :
                       :( CoVectorValue{$Dcov}(($(comps...),)) )
  quote
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_flat(v)
    (; g) = md
    _check_metric_dim(g, $D)
    $result
  end
end

const ♭ = flat

"""
    flat_flat(v::VectorValue{D})
    flat_flat(v::MultiValue{Tuple{D}})
    ♭ᴱ(v)

Musical isomorphism ♭ in Euclidean space, that is [`flat`](@ref) with the
identity metric: the covector (1-form) with the same components as the
contravariant vector `v`. Unlike [`to_1form`](@ref), `v` may already carry its
contravariant variance, and the result always carries the covariant one.
"""
@generated function flat_flat(v::MultiValue{Tuple{Dv},Tv}) where {Dv,Tv}
  Dv isa Con || return :(throw(ArgumentError("flat_flat undefined for 1-tensor that is not contravariant, got $Dv")))
  D = Int(Dv)
  Dcov = inverse_variance(Dv)
  iszero(D) && return :( zero(CoVectorValue{$Dcov,$Tv}) )

  comps = [:( v[$i] ) for i in 1:D]
  :( CoVectorValue{$Dcov}(($(comps...),)) )
end
const ♭ᴱ = flat_flat

"""
    sharp(ω::MultiValue{Tuple{D}}, metric_args...)
    ♯(ω, metric_args...)

Musical isomorphism `♯`: raise the index of the covector (1-form) `ω` with the
inverse metric `g_inv`, resulting into a contravariant vector.

`metric_args` contains the metric information, see [`parse_metric_data`](@ref).
`sharp` uses `g_inv`.

See also [`flat_sharp(ω)`](@ref) for the Euclidean version, accessible by
passing `LinearAlgebra.I`.
"""
@generated function sharp(ω::MultiValue{Tuple{Dw},Tw}, metric_args::Vararg{_MetricArg}) where {Dw,Tw}
  Dw isa Cov || return :(throw(ArgumentError("sharp undefined for 1-tensor that is not covariant, got $Dw")))
  D = Int(Dw)
  Dcon = inverse_variance(Dw)

  # An empty result still promotes with the scalar type of the metric.
  comps = [Expr(:call, :+, (:(g_inv[$i,$j] * indep_comp_getindex(ω,$j)) for j in 1:D)...) for i in 1:D]
  result = iszero(D) ? :( zero(ContraVectorValue{$Dcon,Base.promote_op(*,eltype(g_inv),$Tw)}) ) :
                       :( ContraVectorValue{$Dcon}(($(comps...),)) )
  quote
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_sharp(ω)
    (; g_inv) = md
    _check_metric_dim(g_inv, $D)
    $result
  end
end

const ♯ = sharp

"""
    flat_sharp(ω::MultiValue{Tuple{D}})
    ♯ᴱ(ω)

Musical isomorphism ♯ in Euclidean space, that is [`sharp`](@ref) with the
identity metric: the contravariant vector with the same components as the
covector (1-form) `ω`. Unlike [`from_1form`](@ref), the result carries its
contravariant variance rather than being un-varianced.
"""
@generated function flat_sharp(ω::MultiValue{Tuple{Dw},Tw}) where {Dw,Tw}
  Dw isa Cov || return :(throw(ArgumentError("flat_sharp undefined for 1-tensor that is not covariant, got $Dw")))
  D = Int(Dw)
  Dcon = inverse_variance(Dw)
  iszero(D) && return :( zero(ContraVectorValue{$Dcon,$Tw}) )

  comps = [:( indep_comp_getindex(ω,$i) ) for i in 1:D]
  :( ContraVectorValue{$Dcon}(($(comps...),)) )
end
const ♯ᴱ = flat_sharp


# Application of tensors to vectors and co-vectors (contractions)

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

Pointwise pullback `φ^*ω` of a K-form by a map `φ : ℝᴰᵐ → ℝᴰⁿ` with
transposed Jacobian `Jt = ∇φ`, that is `Jt[j,i] = ∂φⁱ/∂ξʲ`.

The domain index `1≤j≤m` of `Jt` is covariant and its image index `1≤i≤n`
contravariant, so `Jt` has size `(Cov(Dm), Con(Dn))`, the result is an
`ExteriorFormValue{K,Cov(Dm)}`.

A Field-level method for lazy pullbacks is provided in `Gridap.Fields`.
"""
@generated function pullback(ω::ExteriorFormValue{K,Dw,Tw}, Jt::TensorValue{Dm,Dn,Tj}) where {K,Dw,Dm,Dn,Tw,Tj}
  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  can_contract(Dn, Dw) || return :(throw(DimensionMismatch(
    "cannot pullback a form of dimension $Dw by a transposed jacobian of image (second) dimension $Dn, indices cannot contract"
  )))
  T = Base.promote_op(*, Tw, Tj)

  m, n = Int(Dm), Int(Dn)
  (iszero(binomial(n, K)) || iszero(binomial(m, K))) &&
    return :( zero(ExteriorFormValue{$K,$Dm,$T}) )

  c_I = sorted_combinations(n, K)   # K-combins in image (ambient) space
  c_J = sorted_combinations(m, K)   # K-combins in domain (chart) space

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
    koszul(x::MultiValue{Tuple{D}}, ω::ExteriorFormValue)
    koszul(x::VectorValue, ω::ExteriorFormValue)

Koszul differential ``κ_x(ω) = ι_x(ω)``, the interior product with the position
vector `x`, that is either a contravariant vector or a VectorValue.
"""
function koszul(x::MultiValue{Tuple{D}}, ω::ExteriorFormValue{K,Dw}) where {K,D,Dw}
  @assert(D isa Con,
         "the Koszul differential contracts the position vector of a point, whose \
         index is contravariant, got a vector of index $D.")
  interior_product(x, ω)
end

# ============================================================
# Inner product of alternating tensors
#
# For K-forms:
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
    alt_inner(a::ExteriorFormValue{K,D},   b::ExteriorFormValue{K,D},   metric_args...)
    alt_inner(a::ExteriorVectorValue{K,D}, b::ExteriorVectorValue{K,D}, metric_args...)
    ⨟(a, b, metric_args...)

Inner product ``(a|b)`` of two `K`-forms or two `K`-vectors, the scalar such
that, for two `K`-forms,

``a ∧ ⋆ b = √|det(g)| (a|b) dx¹∧…∧dx^D``.

`metric_args` contains the metric information, see [`parse_metric_data`](@ref).
`alt_inner` uses the metric with opposite variance to its arguments: a `K`-form
pair is contracted through `g_inv`, raising the indices of `b` as `hodge_star`
does, and a `K`-vector pair through `g`, lowering them instead.

See also [`flat_alt_inner`](@ref) for the Euclidean version, accessible by
passing `LinearAlgebra.I`.
"""
@generated function alt_inner(a::SkewSymKTensorValue{V,K,Dw,Ta}, b::SkewSymKTensorValue{V,K,Dw,Tb},
                              metric_args::Vararg{_MetricArg}) where {V,K,Dw,Ta,Tb}
  V === Int && return :(@unreachable $("alt_inner requires a K-tensor with variance, got $V index"))
  D = Int(Dw)
  c_K = sorted_combinations(D, K)
  # raising or lowering b's indices as `hodge_star` does, depending on V
  _g_or_ginv = V === Cov ? :g_inv : :g
  # A `K > D` pair has no components; its zero still promotes with the metric.
  T = :( Base.promote_op(*, $Ta, $Tb, eltype($_g_or_ginv)) )
  terms = [:( indep_comp_getindex(a,$m) *
              $(_sum_expr([:( $(_minor_expr(_g_or_ginv, I, J)) * indep_comp_getindex(b,$n) )
                           for (n, J) in enumerate(c_K)], T)) )
           for (m, I) in enumerate(c_K)]
  # (a|b) = Σ_I a_I b^I
  quote
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_alt_inner(a, b)
    (; $_g_or_ginv) = md
    _check_metric_dim($_g_or_ginv, $D)
    $(_sum_expr(terms, T))
  end
end

const ⨟ = alt_inner

"""
    flat_alt_inner(ω::ExteriorFormValue{K,D},   η::ExteriorFormValue{K,D})
    flat_alt_inner(ω::ExteriorVectorValue{K,D}, η::ExteriorVectorValue{K,D})
    ⨟ᴱ(ω,η)
    ω ⨟ᴱ η

Inner product ``(ω|η)`` of two `K`-forms or -vectors in Euclidean space, that is
[`alt_inner`](@ref) with the identity (inverse) metric.

This differs from the tensor inner product by a factor `K!`, see [`flat_inner`](@ref).
"""
@generated function flat_alt_inner(
  a::SkewSymKTensorValue{V,K,D,Ta}, b::SkewSymKTensorValue{V,K,D,Tb}) where {V,K,D,Ta,Tb}

  T = Base.promote_op(*,Ta,Tb)
  L = num_indep_components(a)
  _sum_expr([:(indep_comp_getindex(a,$i) * indep_comp_getindex(b,$i)) for i in 1:L], T)
end
const ⨟ᴱ = flat_alt_inner

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
  factorial(K) * flat_alt_inner(a, b)
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


# ============================================================
# Low level differential operators from partial derivatives
#
# Each of these consumes a tensor of partial derivatives of its arguments: the
# `N` derivation indices come first and are covariant.
# ============================================================

# Exterior derivative

"""
    partials_to_exteriorderivative(∂ω)
    partials_to_exteriorderivative(∂f)

Value of the exterior derivative `dω` of a K-form `ω` at `x`, from its first
order partial derivatives `∂ω = partial_derivatives(ω)(x)`. The second method
takes the partial derivatives (gradient) of a scalar function `f`.

`dω` is the alternation of `∂ω` scaled by `K+1`, so that

    dω[i_0, …, i_K] = ∑_a (-1)^a ∂ω[i_a, i_0, …, î_a, …, i_K]

which is a `ExteriorFormValue{K+1,D}`.
"""
@generated function partials_to_exteriorderivative(
  ∂ω::ProductTensorValue{<:MultiValue{Tuple{Dw}},ExteriorFormValue{K,Dw},T}) where {K,Dw,T}

  D  = Int(Dw)
  Kp1 = K + 1
  ncomp = binomial(D, Kp1)
  iszero(ncomp) &&
    return :( zero(ExteriorFormValue{$Kp1,$Dw,$(Base.promote_op(-,T,T))}) )

  # (dω)_I = Σ_q (-1)^{q-1} ∂_{I[q]} ω_{I∖I[q]}
  comps = Vector{Expr}(undef, ncomp)
  for (w,I) in enumerate(sorted_combinations(D, Kp1))
    terms = Vector{Expr}(undef, Kp1)
    for (q, I_sub_q_id) in enumerate(sub_combinations_ids(I,D))
      sgn = isodd(q) ? 1 : -1
      i = (I_sub_q_id-1)*D + I[q]
      terms[q] = _signed(sgn, :( indep_comp_getindex(∂ω, $i) ))
    end
    comps[w] = _sum_expr(terms, T)
  end

  :( ExteriorFormValue{$Kp1,$Dw}(($(comps...),)) )
end

function partials_to_exteriorderivative(∂f::MultiValue{Tuple{Dw},T,1}) where {Dw,T}
  Dw isa Cov || throw(ArgumentError(
    "the exterior derivative expect covariant (`Cov`) partial derivatives, got a derivation index of variance $(typeof(Dw))"))
  ExteriorFormValue{1,Dw,T}(Tuple(∂f))
end


# Lie derivative

"""
    partials_to_liederivative(u, v, ∂u, ∂v)
    partials_to_liederivative(u, ω::ExteriorFormValue, ∂u, ∂ω)

Lie derivative `𝓛_u v` of a field `v` along the contravariant vector field `u`,
from both fields and their first order partial derivatives. `v` may be a
function, a vector, a covector or a K-form, and `𝓛_u v` has its variance and
degree.

On a K-form, Cartan's magic formula `𝓛_u ω = d(ι_u ω) + ι_u(dω)` expands to

    𝓛_u ω[i_1, …, i_K] = ∑_j u[j]*∂ω[j, i_1, …, i_K]
                       + ∑_a ∑_j ω[i_1, …, j, …, i_K]*∂u[i_a, j]

the second sum replacing `i_a` by `j`. On a vector field it is the Lie bracket
`𝓛_u v = [u, v]`, on a function, 𝓛_u f = u⋅∂f.
"""
function partials_to_liederivative(v::MultiValue, ω::Number, ∂v::MultiValue, ∂ω::MultiValue)
  @notimplemented
end

# Lie derivative of a scalar function
function partials_to_liederivative(
  u::MultiValue{Tuple{Du}}, f::Number,
  ∂u::MultiValue{Tuple{Dp,Du}}, ∂f::MultiValue{Tuple{Dp}}) where {Du,Dp}

  Du isa Con || throw(ArgumentError(
    "Lie derivative requires a contravariant first argument"))
  Dp isa Cov || throw(ArgumentError(
    "partial derivatives should have a covariant index, got $Dp"))

  # 𝓛ᵤf = u(f) = df(u) = uʲ ∂ⱼf
  return u⋅∂f
end

# Lie derivative of a vector field
function partials_to_liederivative(
  u::MultiValue{Tuple{D}}, v::MultiValue{Tuple{D}},
  ∂u::MultiValue{Tuple{Dp,D}}, ∂v::MultiValue{Tuple{Dp,D}}) where {D,Dp}

  D == Dp || throw(DimensionMismatch(
    "the arguments of Lie derivative must all be expressed in the same coordinates, got incompatible dimensions $D ≠ $Dp"))
  D isa Con || throw(ArgumentError(
    "Lie derivative requires a contravariant first argument"))
  Dp isa Cov || throw(ArgumentError(
    "partial derivatives should have a covariant index, got $Dp"))

  # (𝓛ᵤv)ⁱ = [u, v]ⁱ = uʲ ∂ⱼvⁱ - vʲ ∂ⱼuⁱ
  return u⋅∂v - v⋅∂u
end

# Lie derivative of a covector
function partials_to_liederivative(
  u::MultiValue{Tuple{Du}}, ω::MultiValue{Tuple{Dv}},
  ∂u::MultiValue{Tuple{Dv,Du}}, ∂ω::MultiValue{Tuple{Dv,Dv}}) where {Du,Dv}

  Du == Dv || throw(DimensionMismatch(
    "the arguments of Lie derivative must all be expressed in the same coordinates, got incompatible dimensions $Du ≠ $Dv"))
  Du isa Con || throw(ArgumentError(
    "Lie derivative requires a contravariant first argument"))
  Dv isa Cov || throw(ArgumentError(
    "the partial derivatives for the lie derivative should have a covariant index, got $Dv"))

  # (𝓛ᵤω)ᵢ =  uʲ ∂ⱼωᵢ + ∂ᵢuʲ ωⱼ
  return u⋅∂ω + ∂u⋅ω
end

# Lie derivative of a k-form
@generated function partials_to_liederivative(
  u::MultiValue{Tuple{Du},Tu}, ω::ExteriorFormValue{K,Dw,Tw},
  ∂u::MultiValue{Tuple{Dw,Du},Tdu},
  ∂ω::ProductTensorValue{<:MultiValue{Tuple{Dw}},ExteriorFormValue{K,Dw},Tdw}
  ) where {K,Du,Dw,Tu,Tw,Tdu,Tdw}

  Du isa Con || return :(throw(ArgumentError(
    $("the Lie derivative requires a contravariant vector field, got a first argument of size ($Du,)"))))
  Du == Dw || return :(throw(DimensionMismatch(
    $("the vector field and the form must live in the same dimension, got $(Int(Du)) and $(Int(Dw))"))))

  D = Int(Dw)
  T = Base.promote_op(+, Base.promote_op(*,Tu,Tdw), Base.promote_op(*,Tdu,Tw))
  ncomp = binomial(D, K)
  iszero(ncomp) && return :( zero(ExteriorFormValue{$K,$Dw,$T}) )

  comps = Vector{Expr}(undef, ncomp)
  for (n,I) in enumerate(sorted_combinations(D, K))
    # uʲ ∂ⱼω_I, the transport of ω along the flow of u
    terms = Expr[:( u[$j] * ∂ω[$j, $(I...)] ) for j in 1:D]
    # ∂_{i_a}uʲ ω_{i_1…j…i_K}, which vanishes when j already sits in I∖I[a], so
    # for each a the sum runs over the complement of that sub-combination
    for a in 1:K
      Ia = ntuple(i -> I[i + Int(i≥a)], K-1)
      Iacomp = combination_complement(Ia, D)
      for (q, l) in enumerate(sup_combinations_ids(Ia, D))
        j = Iacomp[q]
        sgn = sorting_sign(I[1:a-1]..., j, I[a+1:end]...)
        push!(terms, _signed(sgn, :( ∂u[$(I[a]), $j] * indep_comp_getindex(ω, $l) )))
      end
    end
    comps[n] = _sum_expr(terms, T)
  end

  :( ExteriorFormValue{$K,$Dw}(($(comps...),)) )
end

# TODO later, Lie derivative of arbitrary tensor with no (or unknown) symmetries

# codifferential

"""
    partials_to_codifferential(∂ω)

Value of the flat codifferential `δω` of a K-form `ω` at `x`, a (K-1)-form, from
its first order partial derivatives `∂ω = partial_derivatives(ω)(x)`.

`δω` contracts the derivation index of `∂ω` with the first index of `ω`,

    δω[i_1, …, i_{K-1}] = -∑_j ∂ω[j, j, i_1, …, i_{K-1}]
"""
@generated function partials_to_codifferential(
  ∂ω::ProductTensorValue{<:MultiValue{Tuple{Dw}},ExteriorFormValue{K,Dw},T}) where {K,Dw,T}

  K >= 1 || return :(@unreachable "the codifferential requires K ≥ 1")
  D = Int(Dw)
  Km1 = K - 1
  ncomp = binomial(D, Km1)
  iszero(ncomp) &&
    return :( zero(ExteriorFormValue{$Km1,$Dw,$(Base.promote_op(-,T,T))}) )

  # (δω)_I = -Σ_j ∂_j ω_{j I}. Only j ∉ I contributes, the other terms vanishing
  # by antisymmetry, so the sum runs over the complement of I.
  comps = Vector{Expr}(undef, ncomp)
  for (n,I) in enumerate(sorted_combinations(D, Km1))
    Icomp = combination_complement(I, D)
    terms = Vector{Expr}(undef, length(Icomp))
    for (q, I_sup_q_id) in enumerate(sup_combinations_ids(I, D))
      j = Icomp[q]
      i = (I_sup_q_id-1)*D + j
      terms[q] = _signed(-sorting_sign(j, I...), :( indep_comp_getindex(∂ω, $i) ))
    end
    comps[n] = _sum_expr(terms, T)
  end

  :( ExteriorFormValue{$Km1,$Dw}(($(comps...),)) )
end

# TODO later, codifferential in non-euclidean coordinates
