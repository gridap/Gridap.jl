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
alternation(a::MultiValue{Tuple{}}) = throw(ArgumentError(
  "alternation undefined for $(typeof(a)): it has no dimension nor variance."))
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
# ============================================================

# The step between two components that differ only by index `k`, in the flat,
# column-major layout of a tensor of size `sz`.
_stride(sz, k) = prod(Int.(sz[1:k-1]), init=1)

# Checks that the non-empty index list `ks` of `a` has the variance
# `expected_Dim` and that one metric can act on all of it, and prepares result
# data: the type of the result, the stride of each index of `ks` in the flat
# column-major layout, the dimension they share and the component count.
#
# Returns `(nothing, expr)` with `expr` an `Expr` that the calling `@generated`
# function must return as is to propagate the error.
function _parse_musical_indices(a, ks, expected_Dim, op)
  sz = size(a)
  N = length(sz)
  variance = expected_Dim === Con ? "contravariant" : "covariant"

  all(k -> k isa Int && 1 <= k <= N, ks) || return nothing, :( throw(ArgumentError(
    $("$op undefined for the indices $ks of a tensor of order $N"))) )
  allunique(ks) || return nothing, :( throw(ArgumentError(
    $("$op cannot act twice on one index, got the indices $ks"))) )
  all(k -> sz[k] isa expected_Dim, ks) || return nothing, :( throw(ArgumentError(
    $("$op undefined for the indices $ks of a tensor of size $sz, they are not all $variance"))) )

  # One metric acts on every lowered or raised index, so they all share the
  # dimension it is given for.
  D = Int(sz[first(ks)])
  all(k -> Int(sz[k]) == D, ks) || return nothing, :( throw(DimensionMismatch(
    $("the indices $ks of a tensor of size $sz have different dimensions, so no \
       single metric can act on them all"))) )


  Sr = ntuple(i -> i in ks ? inverse_variance(sz[i]) : sz[i], N)
  (V = _tensor_base_type_expr(N, Sr), ss = map(k -> _stride(sz, k), ks),
   D = D, L = length(a)), nothing
end

# Lowering or raising of the indices of strides `ss`, each of dimension `D`, of
# the `L` components of `a` by the metric named `g_or_ginv`.
function _musical_iterations_expr(V, ss, D, L, g_or_ginv, Ta)
  iterations = Expr(:block)
  src, dst = :x0, :x1
  for s in ss
    push!(iterations.args, :( _data_single_metric_contract!($dst, $src, $g_or_ginv, Val($s), Val($D)) ))
    src, dst = dst, src
  end
  quote
    T = Base.promote_op(*, $Ta, eltype($g_or_ginv))
    # these flattened reusable buffers enable minimizing the data shuffles
    x0 = MVector{$L,T}(Tuple(get_array(a)))
    x1 = MVector{$L,T}(undef)
    $iterations
    $V(Tuple($src))
  end
end

# Lower or raise the `Dk` dimentional index `k` of tensor data `x`
#
#   y[i1,…,ik,…,iN] = Σ_j M[ik,j] * x[i1,…,j,…,iN]
#
# `y` and `x` must be distinct buffers. The stride `s` of `k` is the number of
# components between two that differ only by ±1 their `k`th index.
@inline function _data_single_metric_contract!(
  y::MVector{L,T}, x::StaticVector{L}, g_or_ginv, ::Val{s}, ::Val{Dk}
) where {L,T,s,Dk}

  # TODO: check in marking the whole block @inbounds makes it faster
  for o in 0:(L ÷ (s*Dk) - 1)
    base = o*s*Dk
    for p in 1:s
      q = base + p
      for i in 1:Dk
        acc = zero(T)
        for j in 1:Dk
          acc += g_or_ginv[i,j] * x[q + (j-1)*s]
        end
        y[q + (i-1)*s] = acc
      end
    end
  end
  nothing
end

const _SymValue = Union{AbstractSymMultiValue,SymFourthOrderTensorValue}

# unlike MultiValue(SArray(a)), this preserves variance and return ScalarValue for 0-tensors
@generated function full_rank(a::MultiValue{S,T,N}) where {S,T,N}
  iszero(N) && return :( ScalarValue{$T}((a[1],)) ) # V is :() for O order tensor
  V = _tensor_base_type_expr(N, size(a))
  iszero(length(a)) && return :( zero($V{$T}) )
  :( $V(Tuple(get_array(a))) )
end

# Whether a musical operator flipping the indices `ks`, of variance `V`, of a
# tensor of type `F` returns a tensor with the same symmetries.
function _musical_preserves_symmetry(F, ks, ::Type{V}) where V
  isempty(ks) && return true
  F <: SkewSymKTensorValue && return ks === ntuple(identity, length(size(F)))
  F <: _SymValue || return true

  sz = size(F)
  N = length(sz)
  # `ks` holds every index once: it has `N` of them and each of `1:N` is in it.
  length(ks) === N && all(k -> k isa Int, ks) && all(∈(ks), 1:N) || return false
  # all indices are iddentical
  first(sz) isa V && all(d -> d === first(sz), sz)
end

_musical_symmetric_type(::Type{<:AbstractSymTensorValue{D}}) where D =
  SymTensorValue{inverse_variance(D)}
_musical_symmetric_type(::Type{<:SymNTensorValue{S}}) where S =
  SymNTensorValue{inverse_variance(S)}
_musical_symmetric_type(::Type{<:SymFourthOrderTensorValue{D}}) where D =
  SymFourthOrderTensorValue{inverse_variance(D)}


# `b_I = Σ_J G[i₁,j₁] … G[i_N,j_N] a_J`, the congruence by the symmetric matrix
# `G` that every index of `a` is flipped by.
#
# The `N-1` first indices are contracted by the mode-product ping-pong the
# generic methods use, whose loop LLVM vectorizes. The last one is contracted
# straight into the components the result stores, which costs
# `num_indep_components·D` products where one more pass would cost `D^(N+1)`.
# No component of the result but those is formed and no intermediate
# `MultiValue` is constructed; the buffer the passes work in still holds the
# whole `D^N` of the argument.
#
# Running the last index through the loop as well and only reading the stored
# components off the buffer is the other option, and a slower one everywhere but
# on a `SymFourthOrderTensorValue{con3}`, where the two are within 1%:
# 23.8 against 18.5 ns on a `SymTensorValue{con3}`, 33.6 against 29.1 ns on a
# `SymTensorValue{con4}` and 64.0 against 36.8 ns on a
# `SymNTensorValue{Tuple{con3,con3,con3}}`.
@generated function _musical_congruence_kernel(a::MultiValue{S,Ta,N}, g_or_ginv) where {S,Ta,N}
  D = Int(size(a)[1])
  L = D^N
  Vr = _musical_symmetric_type(a)
  T = :( Base.promote_op(*, eltype(g_or_ginv),Ta) )
  iszero(L) && return :( zero($Vr{$T}) )

  passes = Expr(:block)
  src, dst = :x0, :x1
  for k in 1:N-1
    push!(passes.args,
          :( _data_single_metric_contract!($dst, $src, g_or_ginv, Val($(D^(k-1))), Val($D)) ))
    src, dst = dst, src
  end

  # Position of the component `(I₁, …, I_{N-1}, j)` of the buffer the passes
  # leave the `N-1` first indices flipped in, the last one having stride D^(N-1).
  pos(I, j) = 1 + sum(((I[k]-1)*D^(k-1) for k in 1:N-1), init=0) + (j-1)*D^(N-1)
  comps = [ _sum_expr([ :( g_or_ginv[$(I[N]),$j] * y[$(pos(I,j))] ) for j in 1:D ], :T)
            for I in _stored_comp_cartesian_indices(a) ]

  quote
    x0 = MVector{$L,$T}(Tuple(get_array(a)))
    $(N > 1 ? :( x1 = MVector{$L,$T}(undef) ) : nothing)
    $passes
    y = Tuple($src) # for more performent getindex in comps access
    $Vr(($(comps...),))
  end
end

# - `op` ∈ (flat, sharp, flat_flat, flat_sharp)
# - `V<:Dim` the flipped indices variance
# - `ds` is `drop_sym`, enables returning a tensor with less symetries than the input
#
# `@generated`: the index list and the symmetries are both type level, so which
# of the three branches is taken is too. Leaving the choice to the generated
# code would make the return type the union of the symmetric and the full rank
# results, since `_musical_preserves_symmetry` does not constant fold.
@generated function _symmetry_preserving_musical_op(
  op::F, a::_SymValue, ::Val{l}, ::Type{V}, ::Val{ds}, metric_args::Vararg{_MetricArg}) where {F,l,V,ds}

  ks = tuple(l...)
  isempty(ks) && return :( a )
  size_a = size(a)
  flat_res = :( _musical_symmetric_type(typeof(a))(Tuple(a)) )
  g_or_ginv = V === Con ? :g : :g_inv

  _musical_preserves_symmetry(a, ks, V) && return quote
    (op === flat_flat || op === flat_sharp) && return $flat_res
    md = parse_metric_data(metric_args...)
    isidentity(md) && return $flat_res

    (; $g_or_ginv) = md
    D = Int(size(a)[1])
    _check_metric_dim($g_or_ginv, D)
    _musical_congruence_kernel(a, $g_or_ginv)
  end

  # A list that the operation refuses anyway is left to report its own error.
  N = length(size_a)
  if !ds && all(k -> k isa Int && 1 <= k <= N && size_a[k] isa V, ks)
    return :( throw(ArgumentError(
      "This $op operation does not preserve the symmetries of the tensor type (not
        all indices are flipped, or they do not share one variance). If you really
        mean it, pass `drop_sym=Val(true)`, or convert `a` to a tensor type that
        won't loose symmetry through `$op(a, Val($l), …)`."
    )) )
  end

  quote
    op(full_rank(a), Val($l), metric_args...)
  end
end

"""
    flat(v::VectorValue,             metric_args...)
    flat(v::MultiValue{Tuple{D}},    metric_args...)
    flat(b::ExteriorVectorValue,     metric_args...)
    flat(a::MultiValue, ::Val{k},    metric_args...; drop_sym=Val(false))
    flat(a::MultiValue, ::Val{l...}, metric_args...; drop_sym=Val(false))
    ♭(...)

Musical isomorphism `♭` that lowers:
- the only index of the vector `v`, or
- all the indices of the K-vector `b`, or
- the `k`-th index of an arbitrary tensor `a`, or
- the indices `k ∈ l` of `a` all at once.

    (v♭)ᵢ = gᵢⱼ vʲ
    (b♭)_I = Σ_J g_{I₁,J₁} … g_{I₁,J₁} b^J
    (a♭){i₁,…}_iₖ{…,i_N} = Σ_j g_{iₖ,j} a{i₁,…}^j{…,i_N}

The lowered indices must be contravariant. `metric_args` contains the metric
information, see [`parse_metric_data`](@ref). `flat` uses `g`.

A tensor type containing symmetries only keeps them if every index is lowered
to `Cov`. Otherwise, `flat` errors on such types rather than removing the
symmetries, except if `drop_sym=Val(true)` is passed.

For [`ProductTensorValue`](@ref), the previous rule applies to each product
factor types independently. `drop_sym` applies to every factor.

See also [`flat_flat`](@ref) for the Euclidean version, accessible by passing
`LinearAlgebra.I`.
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

@generated function flat(a::MultiValue{S,Ta}, ::Val{l}, metric_args::Vararg{_MetricArg};
                         drop_sym::Val{ds}=Val(false)) where {S,Ta,l,ds}
  ks = tuple(l...) # puts scalar l into a tuple
  isempty(ks) && return :( a )
  data, err = _parse_musical_indices(a, ks, Con, "flat")
  err === nothing || return err
  (; V, ss, D, L) = data

  if iszero(L)
    result = :( zero($V{Base.promote_op(*,eltype(g),$Ta)}) )
  else
    result = _musical_iterations_expr(V, ss, D, L, :g, Ta)
  end

  quote
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_flat(a, Val($l))
    (; g) = md
    _check_metric_dim(g, $D)
    $result
  end
end

function flat(a::SkewSymKTensorValue{Con,K}, ::Val{l},
                         metric_args::Vararg{_MetricArg};
                         drop_sym::Val{ds}=Val(false)) where {K,l,ds}
  ks = tuple(l...)
  isempty(ks) && return a
  ks === ntuple(identity, K) && return flat(a, metric_args...)
  ds || throw(ArgumentError(
    "This flat operation does not preserve the skew-symmetry (not all indices
      are lowered). If you really mean it, pass `drop_sym=Val(true)`, or convert
      `a` to a tensor type that won't loose symmetry through
      `flat(a, Val($l), metric_args...)`."))
  flat(full_rank(a), Val(l), metric_args...)
end

# A congruence preserves the index symmetries, so only the independent
# components of the result are kept.
function flat(a::_SymValue, ::Val{l}, metric_args::Vararg{_MetricArg};
              drop_sym::Val{ds}=Val(false)) where {l,ds}
  _symmetry_preserving_musical_op(flat, a, Val(l), Con, Val(ds), metric_args...)
end

# optimization for skew-symmetry flat:
#   ωᵢ₁…ᵢₖ = Σ_{j₁<…<j_K} det(g[I,J]) ω^j₁…^j_K
@generated function flat(b::ExteriorVectorValue{K,Db,Tb},
                         metric_args::Vararg{_MetricArg}) where {K,Db,Tb}
  D = Int(Db)
  Dcov = Cov(D)
  c = sorted_combinations(D, K)

  T = :( Base.promote_op(*, $Tb, eltype(g)) )
  result = if isempty(c)
    # An empty result still promotes with the scalar type of the metric.
    :( zero(ExteriorFormValue{$K,$Dcov,$T}) )
  else
    begin
      comps = [ _sum_expr([ :( $(_minor_expr(:g, I, J)) * indep_comp_getindex(b,$l) )
                            for (l, J) in enumerate(c) ], T)
                for I in c ]
      :( ExteriorFormValue{$K,$Dcov}(($(comps...),)) )
    end
  end

  quote
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_flat(b)
    (; g) = md
    _check_metric_dim(g, $D)
    $result
  end
end

const ♭ = flat

"""
    flat_flat(v::VectorValue{D})
    flat_flat(v::MultiValue{Tuple{D}})
    flat_flat(b::ExteriorVectorValue,)
    flat_flat(a::MultiValue, ::Val{k};    drop_sym=Val(false))
    flat_flat(a::MultiValue, ::Val{l...}; drop_sym=Val(false))
    ♭ᴱ(...)

Musical isomorphism `♭` in Euclidean space, that is [`flat`](@ref) with the
identity metric

    (♭ᴱ(v))ᵢ = vⁱ
    (b♭ᴱ)_I = b^I
    (a♭ᴱ){i₁,…}_iₖ{…,i_N} = a{i₁,…}^iₖ{…,i_N}
"""
@generated function flat_flat(v::MultiValue{Tuple{Dv},Tv}) where {Dv,Tv}
  Dv isa Con || return :(throw(ArgumentError("flat_flat undefined for 1-tensor that is not contravariant, got $Dv")))
  D = Int(Dv)
  Dcov = inverse_variance(Dv)
  iszero(D) && return :( zero(CoVectorValue{$Dcov,$Tv}) )

  comps = [:( v[$i] ) for i in 1:D]
  :( CoVectorValue{$Dcov}(($(comps...),)) )
end
@generated function flat_flat(a::MultiValue{S,Ta}, ::Val{l};
                              drop_sym::Val{ds}=Val(false)) where {S,Ta,l,ds}
  ks = tuple(l...)
  isempty(ks) && return :( a )
  data, err = _parse_musical_indices(a, ks, Con, "flat_flat")
  err === nothing || return err
  (; V, L) = data
  iszero(L) && return :( zero($V{$Ta}) )
  :( $V(Tuple(get_array(a))) )
end

function flat_flat(a::SkewSymKTensorValue{Con,K}, ::Val{l};
                   drop_sym::Val{ds}=Val(false)) where {K,l,ds}
  ks = tuple(l...)
  isempty(ks) && return a
  ks === ntuple(identity, Val(K)) && return flat_flat(a)
  ds || throw(ArgumentError(
    "flat_flat of a skew-symmetric tensor acts on all of its $K indices at once, \
      so the only index list it accepts is $(ntuple(identity, Val(K))), got $l. \
      Pass `drop_sym=Val(true)` for the full rank result."))
  flat_flat(full_rank(a), Val(l))
end

function flat_flat(a::_SymValue, ::Val{l}; drop_sym::Val{ds}=Val(false)) where {l,ds}
  _symmetry_preserving_musical_op(flat_flat, a, Val(l), Con, Val(ds))
end

# Lowering every index of a skew-symmetric tensor Euclideanly relabels its
# variance, the independent components are those of `b` in the same order.
function flat_flat(b::ExteriorVectorValue{K,Db,Tb}) where {K,Db,Tb}
  Dcov = inverse_variance(Db)
  K > Db && return zero(ExteriorFormValue{K,Dcov,Tb})
  ExteriorFormValue{K,Dcov}(Tuple(b))
end
const ♭ᴱ = flat_flat

"""
    sharp(ω::MultiValue{Tuple{D}},    metric_args...)
    sharp(η::ExteriorFormValue,       metric_args...)
    sharp(a::MultiValue, ::Val{k},    metric_args...; drop_sym=Val(false))
    sharp(a::MultiValue, ::Val{l...}, metric_args...; drop_sym=Val(false))
    ♯(...)

Musical isomorphism `♯` that raises:
- the only index of the covector `ω`, or
- all the indices of the K-form `η`, or
- the `k`-th index of an arbitrary tensor `a`, or
- the indices `k ∈ l` of `a` all at once.

    (ω♯)ⁱ = g_invⁱʲ ωⱼ
    (η♯)^I = g_inv^{IJ} η_J
    (a♯){i₁,…}^iₖ{…,i_N} = Σ_j g_inv^{iₖ,j} a{i₁,…}_j{…,i_N}

The lowered indices must be covariant. `metric_args` contains the metric
information, see [`parse_metric_data`](@ref). `sharp` uses `g_inv`.

A tensor type containing symmetries only keeps them if every index is raised
to `Con`. Otherwise, `sharp` errors on such types rather than removing the
symmetries, except if `drop_sym=Val(true)` is passed.

For [`ProductTensorValue`](@ref), the previous rule applies to each product
factor types independently. `drop_sym` applies to every factor.

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

@generated function sharp(a::MultiValue{S,Ta}, ::Val{l}, metric_args::Vararg{_MetricArg};
                          drop_sym::Val{ds}=Val(false)) where {S,Ta,l,ds}
  ks = tuple(l...)
  isempty(ks) && return :( a )
  data, err = _parse_musical_indices(a, ks, Cov, "sharp")
  err === nothing || return err
  (; V, ss, D, L) = data

  result = if iszero(L)
    :( zero($V{Base.promote_op(*,eltype(g_inv),$Ta)}) )
  else
    _musical_iterations_expr(V, ss, D, L, :g_inv, Ta)
  end

  quote
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_sharp(a, Val($l))
    (; g_inv) = md
    _check_metric_dim(g_inv, $D)
    $result
  end
end

function sharp(a::SkewSymKTensorValue{Cov,K}, ::Val{l},
                          metric_args::Vararg{_MetricArg};
                          drop_sym::Val{ds}=Val(false)) where {K,l,ds}
  ks = tuple(l...)
  isempty(ks) && return a
  ks === ntuple(identity, K) && return sharp(a, metric_args...)
  ds || throw(ArgumentError(
    "This sharp operation does not preserve the skew-symmetry (not all indices
      are raised). If you really mean it, pass `drop_sym=Val(true)`, or convert
      `a` to a tensor type that won't loose symmetry through
      `sharp(a, Val($l), metric_args...)`."))
  sharp(full_rank(a), Val(l), metric_args...)
end

# A congruence preserves the index symmetries, so only the independent
# components of the result are kept.
function sharp(a::_SymValue, ::Val{l}, metric_args::Vararg{_MetricArg};
               drop_sym::Val{ds}=Val(false)) where {l,ds}
  _symmetry_preserving_musical_op(sharp, a, Val(l), Cov, Val(ds), metric_args...)
end

# optimization for skew-symmetric sharp
#   η^i₁…^i_K = Σ_{j₁<…<j_K} det(g_inv[I,J]) η_j₁…_j_K
@generated function sharp(η::ExteriorFormValue{K,Dη,Tη},
                          metric_args::Vararg{_MetricArg}) where {K,Dη,Tη}
  D = Int(Dη)
  Dcon = inverse_variance(Dη)
  c = sorted_combinations(D, K)

  T = :( Base.promote_op(*, $Tη, eltype(g_inv)) )
  result = if isempty(c)
    :( zero(ExteriorVectorValue{$K,$Dcon,$T}) )
  else
    begin
      comps = [ _sum_expr([ :( $(_minor_expr(:g_inv, I, J)) * indep_comp_getindex(η,$l) )
                            for (l, J) in enumerate(c) ], T)
                for I in c ]
      :( ExteriorVectorValue{$K,$Dcon}(($(comps...),)) )
    end
  end

  quote
    md = parse_metric_data(metric_args...)
    isidentity(md) && return flat_sharp(η)
    (; g_inv) = md
    _check_metric_dim(g_inv, $D)
    $result
  end
end

const ♯ = sharp

"""
    flat_sharp(ω::MultiValue{Tuple{D}})
    flat_sharp(η::ExteriorFormValue)
    flat_sharp(a::MultiValue, ::Val{k};    drop_sym=Val(false))
    flat_sharp(a::MultiValue, ::Val{l...}; drop_sym=Val(false))
    ♯ᴱ(...)

Musical isomorphism `♯` in Euclidean space, that is [`sharp`](@ref) with the
identity metric:

    (ω♯ᴱ)ⁱ = ωᵢ
    (η♯ᴱ)^I = η_I
    (a♯ᴱ){i₁,…}^iₖ{…,i_N} = a{i₁,…}_iₖ{…,i_N}
"""
@generated function flat_sharp(ω::MultiValue{Tuple{Dw},Tw}) where {Dw,Tw}
  Dw isa Cov || return :(throw(ArgumentError("flat_sharp undefined for 1-tensor that is not covariant, got $Dw")))
  D = Int(Dw)
  Dcon = inverse_variance(Dw)
  iszero(D) && return :( zero(ContraVectorValue{$Dcon,$Tw}) )

  comps = [:( indep_comp_getindex(ω,$i) ) for i in 1:D]
  :( ContraVectorValue{$Dcon}(($(comps...),)) )
end
@generated function flat_sharp(a::MultiValue{S,Ta}, ::Val{l};
                               drop_sym::Val{ds}=Val(false)) where {S,Ta,l,ds}
  ks = tuple(l...)
  isempty(ks) && return :( a )
  data, err = _parse_musical_indices(a, ks, Cov, "flat_sharp")
  err === nothing || return err
  (; V, L) = data
  iszero(L) && return :( zero($V{$Ta}) )
  :( $V(Tuple(get_array(a))) )
end

function flat_sharp(a::SkewSymKTensorValue{Cov,K}, ::Val{l};
                    drop_sym::Val{ds}=Val(false)) where {K,l,ds}
  ks = tuple(l...)
  isempty(ks) && return a
  ks === ntuple(identity, K) && return flat_sharp(a)
  ds || throw(ArgumentError(
    "flat_sharp of a skew-symmetric tensor acts on all of its $K indices at \
       once, so the only index list it accepts is $(ntuple(identity, K)), got \
       $l. Pass `drop_sym=Val(true)` for the full rank result."))
  flat_sharp(full_rank(a), Val(l))
end

function flat_sharp(a::_SymValue, ::Val{l}; drop_sym::Val{ds}=Val(false)) where {l,ds}
  _symmetry_preserving_musical_op(flat_sharp, a, Val(l), Cov, Val(ds))
end

function flat_sharp(η::ExteriorFormValue{K,Dη,Tη}) where {K,Dη,Tη}
  Dcon = inverse_variance(Dη)
  K > Dη && return zero(ExteriorVectorValue{K,Dcon,Tη})
  ExteriorVectorValue{K,Dcon}(Tuple(η))
end
const ♯ᴱ = flat_sharp


# ============================================================
# Musical isomorphisms on a product tensor
# ============================================================

# Musicals apply on each factor independently:
#
#   ♭(a ⊗ b, (ks_a, ks_b)) = ♭(a, ks_a) ⊗ ♭(b, ks_b)
#
# The symmetry preserving rules are thus applied for each factors indices independently.

# A product stores the independent components of its first factor fastest, so
# `b_rows[lar]` holds, for each `lbr`, the component `(lar, lbr)` of the result.
@generated function _ProductTV_from_second_factor_rows(
  ::Type{A}, b_rows::NTuple{iLa,B}) where {A,B,iLa}

  iszero(iLa) && return :(@unreachable "Non-empty first factor is assumed here")
  comps = [ :( indep_comp_getindex(b_rows[$la], $lb) )
            for lb in 1:num_indep_components(B) for la in 1:iLa ]
  :( ProductTensorValue{$A,$B}(($(comps...),)) )
end

@generated function _ProductTV_musical_apply_second_factor(
  op::F, acols::NTuple{iLb,A}, ::Type{B}, ::Val{ksb},
  metric_args::Vararg{_MetricArg}) where {F,A,iLb,B,ksb}

  iszero(iLb) && return :(@unreachable "Non-empty second factor is assumed here")
  iLa = num_indep_components(A)
  rows = [ :( op($B(($(( :( indep_comp_getindex(acols[$lb], $la) ) for lb in 1:iLb)...),)),
                Val($ksb), metric_args...; drop_sym=Val(true)) )
           for la in 1:iLa ]
  :( _ProductTV_from_second_factor_rows($A, ($(rows...),)) )
end

#
@generated function _ProductTV_musical(
  op::F, p::ProductTensorValue{A,B,Tp}, ::Val{l}, ::Type{V}, ::Val{ds},
  metric_args::Vararg{_MetricArg}) where {F,A,B,Tp,l,V,ds}

  ks = tuple(l...)
  isempty(ks) && return :( p )
  op_name = string(nameof(F.instance))
  _, err = _parse_musical_indices(p, ks, V, op_name)
  err === nothing || return err

  Na = length(size(A))
  # flipped indices of each factor
  ksa, ksb = (filter(<=(Na), ks), map(k -> k - Na, filter(>(Na), ks)))

  if !ds
    wrong_factor = !_musical_preserves_symmetry(A, ksa, V) ? (A, "first")  :
                   !_musical_preserves_symmetry(B, ksb, V) ? (B, "second") : nothing

    if !isnothing(wrong_factor)
      factor, which = wrong_factor
      _msg = "This $op_name operation does not preserve the symmetries of the $which factor
            ($factor) of the product tensor (not all of its indices are flipped, or
             they do not share one variance). If you really mean it, pass
             `drop_sym=Val(true)` to flip the indices of the full rank factor instead, or
             convert that factor to a tensor type that won't loose symmetry through
             `$op_name(a, Val($l), …)`."
      return :( throw(ArgumentError($_msg)) )
    end
  end

  iLa, iLb = num_indep_components(A), num_indep_components(B)
  # No component to flip: the result is the product of the flipped zero factors.
  iszero(iLa*iLb) && return :( product_tensor(
    op(zero($(make_concretetype(A{Tp}))), Val($ksa), metric_args...; drop_sym=Val(true)),
    op(zero($(make_concretetype(B{Tp}))), Val($ksb), metric_args...; drop_sym=Val(true))) )

  # to parse (maybe compute) the correct metric tensor once for all
  _metric_args_expr = if op === typeof(flat)
    :( (; g)     = parse_metric_data(metric_args...); _metric_args = (; g) )
  elseif op === typeof(sharp)
    :( (; g_inv) = parse_metric_data(metric_args...); _metric_args = (; g_inv)  )
  else # flat_*
    :( _metric_args = () )
  end

  # Apply op of the indices of the first factor A
  cols = [ :( op($A(($(( :( ic[$((lb-1)*iLa + la)] ) for la in 1:iLa)...),)),
                Val($ksa), _metric_args...; drop_sym=Val(true)) )
           for lb in 1:iLb ]

  quote
    ic = get_indep_components(p)
    $_metric_args_expr
    _ProductTV_musical_apply_second_factor(op, ($(cols...),), $B, Val($ksb), _metric_args...)
  end
end

function flat(p::ProductTensorValue, ::Val{l}, metric_args::Vararg{_MetricArg};
              drop_sym::Val{ds}=Val(false)) where {l,ds}
  _ProductTV_musical(flat, p, Val(l), Con, Val(ds), metric_args...)
end

function sharp(p::ProductTensorValue, ::Val{l}, metric_args::Vararg{_MetricArg};
               drop_sym::Val{ds}=Val(false)) where {l,ds}
  _ProductTV_musical(sharp, p, Val(l), Cov, Val(ds), metric_args...)
end

function flat_flat(p::ProductTensorValue, ::Val{l};
                   drop_sym::Val{ds}=Val(false)) where {l,ds}
  _ProductTV_musical(flat_flat, p, Val(l), Con, Val(ds))
end

function flat_sharp(p::ProductTensorValue, ::Val{l};
                    drop_sym::Val{ds}=Val(false)) where {l,ds}
  _ProductTV_musical(flat_sharp, p, Val(l), Cov, Val(ds))
end

# The only index of an order 1 product is the only index of its order 1 factor.
flat(p::ProductTensorValue{A,B,T,Tuple{D}}, metric_args::Vararg{_MetricArg}) where {A,B,T,D} =
  _ProductTV_musical(flat, p, Val((1,)), Con, Val(false), metric_args...)

sharp(p::ProductTensorValue{A,B,T,Tuple{D}}, metric_args::Vararg{_MetricArg}) where {A,B,T,D} =
  _ProductTV_musical(sharp, p, Val((1,)), Cov, Val(false), metric_args...)

flat_flat(p::ProductTensorValue{A,B,T,Tuple{D}}) where {A,B,T,D} =
  _ProductTV_musical(flat_flat, p, Val((1,)), Con, Val(false))

flat_sharp(p::ProductTensorValue{A,B,T,Tuple{D}}) where {A,B,T,D} =
  _ProductTV_musical(flat_sharp, p, Val((1,)), Cov, Val(false))


###################################################################
# Application of tensors to vectors and co-vectors (contractions) #
###################################################################

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

The domain index `1≤j≤Dm` of `Jt` is covariant and its image index `1≤i≤Dn`
contravariant, so `Jt` has size `(Cov(Dm), Con(Dn))`, the result is an
`ExteriorFormValue{K,Cov(Dm)}`.
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

"""
    pushforward(v::ExteriorVectorValue{K,Dm}, Jt::TensorValue{Dm,Dn})
    pushforward(v::ExteriorVectorValue{K,Dm}, Jt::TensorValue{Dm,Dn}, iJt::TensorValue)

Pointwise pushforward `φ_*v` of a K-vector by a map `φ : ℝᴰᵐ → ℝᴰⁿ` with
transposed Jacobian `Jt = ∇φ`, that is `Jt[j,i] = ∂φⁱ/∂ξʲ`.

The domain index `1≤j≤Dm` of `Jt` is covariant and its image index `1≤i≤Dn`
contravariant, so `Jt` has size `(Cov(Dm), Con(Dn))`, the result is an
`ExteriorVectorValue{K,Con(Dn)}`.

Unlike [`pullback`](@ref), which contracts `Jt`'s image (second) index against
the form, `pushforward` contracts its domain (first) index against the
K-vector, so it is defined directly from `Jt` with no (pseudo-)inversion. The
three argument form ignores `iJt`, it is accepted so that a caller holding a
precomputed `iJt = pinvJt(Jt)` transports every kind of tensor with one call.
"""
@generated function pushforward(v::ExteriorVectorValue{K,Dw,Tw}, Jt::TensorValue{Dm,Dn,Tj}) where {K,Dw,Dm,Dn,Tw,Tj}
  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  can_contract(Dm, Dw) || return :(throw(DimensionMismatch(
    "cannot pushforward a vector of dimension $Dw by a transposed jacobian of domain (first) dimension $Dm, indices cannot contract"
  )))
  T = Base.promote_op(*, Tw, Tj)

  m, n = Int(Dm), Int(Dn)
  (iszero(binomial(n, K)) || iszero(binomial(m, K))) &&
    return :( zero(ExteriorVectorValue{$K,$Dn,$T}) )

  c_I = sorted_combinations(n, K)   # K-combins in image (ambient) space
  c_J = sorted_combinations(m, K)   # K-combins in domain (chart) space

  # det([∂φ^{I[b]}/∂ξ^{J[a]}]) is the minor of Jt on rows J and columns I
  comps = [_sum_expr([:( $(_minor_expr(:Jt, J, I)) * indep_comp_getindex(v,$l) )
                      for (l, J) in enumerate(c_J)], T)
           for I in c_I]

  :( ExteriorVectorValue{$K,$Dn}(($(comps...),)) )
end

# ============================================================
# Transport of a tensor of type (r,s,0) by a map
#
# Given φ : ℝᴰᵐ → ℝᴰⁿ with transposed Jacobian Jt = ∇φ and iJt = pinvJt(Jt),
# a covariant index is pulled back by Jt and pushed forward by iJt, and a
# contravariant index is pushed forward by Jt and pulled back by iJt:
#
#   (φ^*a)_J = Σ_I a_I · M₁(J₁,I₁) ⋯ M_K(J_K,I_K)
#
# In both directions the matrix of a covariant index is read [out,in] and that
# of a contravariant one [in,out]; only which of Jt and iJt plays each role
# changes, which is what the `cov_mat`/`con_mat` arguments below select.
# ============================================================

# Number of scalar `K`-factor product terms in the direct, fully unrolled
# expansion of an order-`K` transport between dimensions `n_in` and `n_out`:
# `n_out^K` output components, each a sum of `n_in^K` products.
_pullpush_num_terms(n_out::Int, n_in::Int, K::Int) = (n_out*n_in)^K

# Above this many terms, the direct expansion makes its own @generated compile
# time impractical: that time grows much faster than linearly in the term count
# (729 terms: 0.2s, 1296: 1.8s, 2401: 4.5s, 4096: 10.3s). The threshold trades
# a bit of runtime speed -- the direct expansion is ~3x faster than the
# iterative one -- for keeping the one-time compile cost of any new tensor
# shape under ~2s.
const _direct_pullpush_term_threshold = 2000

# `iJt = pinvJt(Jt)` is a true inverse only for a square `Jt`; for
# `Int(Dm) < Int(Dn)` it is a left inverse, exact on an object carried by the
# image tangent space, and `Int(Dm) > Int(Dn)` has no such inverse. Returns
# `nothing` when `Jt` admits the transport, else the `Expr` the calling
# `@generated` function returns to raise at runtime.
function _inverse_pullpush_error(Dm, Dn, op)
  Int(Dm) <= Int(Dn) && return nothing
  msg = "$op by the inverse map requires Int(Dm) ≤ Int(Dn) (got Dm=$Dm,
    Dn=$Dn), see `pinvJt`"
  :(throw(ArgumentError($msg)))
end

# `pinvJt` swaps the position and the variance of both indices, so this is the
# only shape an `iJt` of a `Jt` of size `(Dm, Dn)` can have. The comparison is
# by identity: `==` between two `Dim`s compares their value and ignores their
# variance, so it would accept a transposed `iJt`.
function _pinvJt_shape_error(Dm, Dn, Di1, Di2)
  (Di1 === inverse_variance(Dn) && Di2 === inverse_variance(Dm)) && return nothing
  msg = "iJt must equal pinvJt(Jt): expected size
    ($(inverse_variance(Dn)), $(inverse_variance(Dm))), got ($Di1, $Di2)"
  :(throw(ArgumentError($msg)))
end

# Validates that every index of `a` (of size `sz`) can be transported, and
# classifies each as covariant or contravariant. `Dc` is the contravariant
# dimension that a covariant index of `a` must contract with: the image
# dimension of `Jt` for a pullback, the flipped domain dimension for a
# pushforward. Returns `(is_cov, nothing)` on success, or `(nothing, expr)`
# with `expr` an `Expr` the calling `@generated` function should return
# immediately to propagate the error.
function _pullpush_validate_indices(sz, Dc, op)
  is_cov = Bool[]
  for k in eachindex(sz)
    Sk = sz[k]
    if Sk isa Int
      msg = "$op requires every index of `a` to be covariant or contravariant,
        got index $k of dimension $Sk"
      return nothing, :(throw(ArgumentError($msg)))
    else # Sk is varianced
      if !can_contract(Sk isa Cov ? Dc : inverse_variance(Dc), Sk)
        msg = "cannot $op a tensor whose index $k has dimension $Sk, every
          index must be covariant or contravariant of dimension $(Int(Dc))"
        return nothing, :(throw(DimensionMismatch($msg)))
      end
    end
    push!(is_cov, Sk isa Cov)
  end
  is_cov, nothing
end

_pullpush_eltype(Ta, is_cov, Tcov, Tcon) =
  all(is_cov) ? Base.promote_op(*, Ta, Tcov)        :
  any(is_cov) ? Base.promote_op(*, Ta, Tcov, Tcon)  :
                Base.promote_op(*, Ta, Tcon)

# Fully unrolled transport body, with `a` and the result dense (not
# alternating), so every input index tuple `I ∈ (1:n_in)^K` and every output
# tuple `O ∈ (1:n_out)^K` contributes. Assumes `a` and the matrices named by
# `cov_mat` and `con_mat` are bound in the generated code's scope.
function _pullpush_direct_body(is_cov, n_out, n_in, K, T, V, cov_mat, con_mat)
  o_ci = CartesianIndices(ntuple(_ -> n_out, K))
  i_ci = CartesianIndices(ntuple(_ -> n_in, K))
  comps = [_sum_expr(
             [Expr(:call, :*, :( a[$I] ),
                   (is_cov[k] ? :( $cov_mat[$(O[k]),$(I[k])] ) :
                                :( $con_mat[$(I[k]),$(O[k])] ) for k in 1:K)...)
              for I in i_ci],
             T)
           for O in o_ci]
  :( $V($(comps...)) )
end

# Same transport, computed by contracting one index of `a` at a time with
# `tensor_contraction` instead of unrolling the full `K`-way sum, to keep the
# generated code size (and so compile time) bounded for large `K`.
# `tensor_contraction` always appends the newly contracted index at the end of
# its result, so repeatedly contracting the *current* tensor's first index,
# processing `a`'s original indices 1,…,K in order, leaves the final indices
# in that same original order.
function _pullpush_iterative_body(is_cov, K, cov_mat, con_mat)
  body = Expr(:block, :( c = a ))
  for k in 1:K
    M, slot = is_cov[k] ? (cov_mat, 2) : (con_mat, 1)
    push!(body.args, :( c = tensor_contraction(c, $M, Val(1), Val($slot)) ))
  end
  push!(body.args, :c)
  body
end

_pullpush_core_body(is_cov, n_out, n_in, K, T, V, cov_mat, con_mat) =
  _pullpush_num_terms(n_out, n_in, K) > _direct_pullpush_term_threshold ?
    _pullpush_iterative_body(is_cov, K, cov_mat, con_mat) :
    _pullpush_direct_body(is_cov, n_out, n_in, K, T, V, cov_mat, con_mat)

"""
    pullback(a::MultiValue{S}, Jt::TensorValue{Dm,Dn})
    pullback(a::MultiValue{S}, Jt::TensorValue{Dm,Dn}, iJt::TensorValue)

Pointwise pullback `φ^*a` of a tensor `a` of type `(r,s)` (`r` covariant and
`s` contravariant indices, in the order given by `S`) by a map `φ : ℝᴰᵐ → ℝᴰⁿ`
with transposed Jacobian `Jt = ∇φ`, that is `Jt[j,i] = ∂φⁱ/∂ξʲ`:

    (φ^*a)_{J} = Σ_{I} a_{I} · M₁(J₁,I₁) ⋯ M_K(J_K,I_K)

    M_k(j,i) = Jt[j,i]    if the `k`-th index of `a` is covariant
               iJt[i,j]   if the `k`-th index of `a` is contravariant

with `iJt = pinvJt(Jt)`. The domain index `1≤j≤Dm` of `Jt` is covariant and
its image index `1≤i≤Dn` contravariant, so `Jt` has size
`(Cov(Dm), Con(Dn))`. The result is a tensor of the same `(r,s)` type, with every
index now of dimension `Dm`.

Unlike [`pullback`](@ref)`(::ExteriorFormValue, ::TensorValue)`, `a` need not
be alternating: no symmetry of `a` is assumed or preserved, and the result is
a dense `CoVectorValue`/`ContraVectorValue`, `TensorValue` or `NTensorValue`
(for `K = 1`, `2` or `≥ 3` respectively).

A contravariant index is only ever pulled back exactly when `Int(Dm)` equals
`Int(Dn)`: for `Int(Dm) < Int(Dn)` (e.g. `φ` the embedding of a lower
dimensional domain), `iJt` is a pseudo-inverse and a contravariant index gets
projected onto the image tangent space instead, see [`pinvJt`](@ref).
`Int(Dm) > Int(Dn)` is not supported for a contravariant index.

Above `_direct_pullpush_term_threshold` terms, `a` is pulled back index by
index with [`tensor_contraction`](@ref) instead of a single fully unrolled
expression, to keep `@generated` compile times bounded.

If `a` has a contravariant index, `pullback(a,Jt)` computes `iJt = pinvJt(Jt)`
itself on every call. When pulling back many tensors by the same `Jt` (e.g.
across the quadrature points of one cell), pass a precomputed
`iJt = pinvJt(Jt)` explicitly with the 3-argument form to avoid recomputing it
each time; `iJt` is trusted to equal `pinvJt(Jt)` and is not checked against
it, only against its expected size.
"""
@generated function pullback(a::MultiValue{S,Ta}, Jt::TensorValue{Dm,Dn,Tj}) where {S,Ta,Dm,Dn,Tj}
  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  sz = size(a)
  isempty(sz) && return :( a )
  is_cov, err = _pullpush_validate_indices(sz, Dn, "pullback")
  err !== nothing && return err

  if !all(is_cov)
    err = _inverse_pullpush_error(Dm, Dn, "pullback")
    err === nothing || return err
    return :( pullback(a, Jt, pinvJt(Jt)) )
  end

  K = length(sz)
  m, n = Int(Dm), Int(Dn)
  T = _pullpush_eltype(Ta, is_cov, Tj, Tj)
  V = _tensor_base_type_expr(K, ntuple(_ -> Dm, K))

  (iszero(m) || iszero(n)) && return :( zero($V{$T}) )

  _pullpush_core_body(is_cov, m, n, K, T, V, :Jt, :iJt)
end

@generated function pullback(
  a::MultiValue{S,Ta}, Jt::TensorValue{Dm,Dn,Tj}, iJt::TensorValue{Di1,Di2,Ti}
  ) where {S,Ta,Dm,Dn,Tj,Di1,Di2,Ti}

  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  err = _pinvJt_shape_error(Dm, Dn, Di1, Di2)
  err === nothing || return err
  sz = size(a)
  isempty(sz) && return :( a )
  is_cov, err = _pullpush_validate_indices(sz, Dn, "pullback")
  err !== nothing && return err

  if !all(is_cov)
    err = _inverse_pullpush_error(Dm, Dn, "pullback")
    err === nothing || return err
  end

  K = length(sz)
  m, n = Int(Dm), Int(Dn)
  T = _pullpush_eltype(Ta, is_cov, Tj, Ti)
  Sr = ntuple(k -> is_cov[k] ? Dm : inverse_variance(Dm), K)
  V = _tensor_base_type_expr(K, Sr)

  (iszero(m) || iszero(n)) && return :( zero($V{$T}) )

  _pullpush_core_body(is_cov, m, n, K, T, V, :Jt, :iJt)
end

"""
    pushforward(a::MultiValue{S}, Jt::TensorValue{Dm,Dn})
    pushforward(a::MultiValue{S}, Jt::TensorValue{Dm,Dn}, iJt::TensorValue)

Pointwise pushforward `φ_*a` of a tensor `a` of type `(r,s,0)` by a map
`φ : ℝᴰᵐ → ℝᴰⁿ` with transposed Jacobian `Jt = ∇φ`:

    (φ_*a)_{I} = Σ_{J} a_{J} · M₁(I₁,J₁) ⋯ M_K(I_K,J_K)

    M_k(i,j) = iJt[i,j]   if the `k`-th index of `a` is covariant
               Jt[j,i]    if the `k`-th index of `a` is contravariant

with `iJt = pinvJt(Jt)`. Every index of `a` must be covariant or contravariant
of dimension `Dm`; the result is a tensor of the same `(r,s,0)` type, with
every index now of dimension `Dn`.

This is the exact mirror of [`pullback`](@ref)`(::MultiValue, ::TensorValue)`,
including the densification of an alternating or symmetric `a`: it is the
contravariant index that needs no (pseudo-)inversion here, while a covariant
index is transported by `iJt` and so requires `Int(Dm) ≤ Int(Dn)`.
"""
@generated function pushforward(a::MultiValue{S,Ta}, Jt::TensorValue{Dm,Dn,Tj}) where {S,Ta,Dm,Dn,Tj}
  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  sz = size(a)
  isempty(sz) && return :( a )
  is_cov, err = _pullpush_validate_indices(sz, inverse_variance(Dm), "pushforward")
  err !== nothing && return err

  if any(is_cov)
    err = _inverse_pullpush_error(Dm, Dn, "pushforward")
    err === nothing || return err
    return :( pushforward(a, Jt, pinvJt(Jt)) )
  end

  K = length(sz)
  m, n = Int(Dm), Int(Dn)
  T = _pullpush_eltype(Ta, is_cov, Tj, Tj)
  V = _tensor_base_type_expr(K, ntuple(_ -> Dn, K))

  (iszero(m) || iszero(n)) && return :( zero($V{$T}) )

  _pullpush_core_body(is_cov, n, m, K, T, V, :iJt, :Jt)
end

@generated function pushforward(
  a::MultiValue{S,Ta}, Jt::TensorValue{Dm,Dn,Tj}, iJt::TensorValue{Di1,Di2,Ti}
  ) where {S,Ta,Dm,Dn,Tj,Di1,Di2,Ti}

  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  err = _pinvJt_shape_error(Dm, Dn, Di1, Di2)
  err === nothing || return err
  sz = size(a)
  isempty(sz) && return :( a )
  is_cov, err = _pullpush_validate_indices(sz, inverse_variance(Dm), "pushforward")
  err !== nothing && return err

  if any(is_cov)
    err = _inverse_pullpush_error(Dm, Dn, "pushforward")
    err === nothing || return err
  end

  K = length(sz)
  m, n = Int(Dm), Int(Dn)
  T = _pullpush_eltype(Ta, is_cov, Ti, Tj)
  Sr = ntuple(k -> is_cov[k] ? inverse_variance(Dn) : Dn, K)
  V = _tensor_base_type_expr(K, Sr)

  (iszero(m) || iszero(n)) && return :( zero($V{$T}) )

  _pullpush_core_body(is_cov, n, m, K, T, V, :iJt, :Jt)
end

# ============================================================
# Alternating transports
#
# A pullback of a K-form and a pushforward of a K-vector contract a minor of
# Jt, the two other transports a minor of iJt = pinvJt(Jt). All four keep the
# alternating type, which the dense (r,s,0) methods above would drop since
# ExteriorFormValue and ExteriorVectorValue are MultiValues too.
# ============================================================

@generated function pushforward(
  v::ExteriorVectorValue, Jt::TensorValue{Dm,Dn,Tj}, iJt::TensorValue{Di1,Di2,Ti}
  ) where {Dm,Dn,Tj,Di1,Di2,Ti}

  err = _pinvJt_shape_error(Dm, Dn, Di1, Di2)
  err === nothing || return err
  :( pushforward(v, Jt) )
end

@generated function pullback(
  ω::ExteriorFormValue, Jt::TensorValue{Dm,Dn,Tj}, iJt::TensorValue{Di1,Di2,Ti}
  ) where {Dm,Dn,Tj,Di1,Di2,Ti}

  err = _pinvJt_shape_error(Dm, Dn, Di1, Di2)
  err === nothing || return err
  :( pullback(ω, Jt) )
end

"""
    pushforward(ω::ExteriorFormValue{K,Dm}, Jt::TensorValue{Dm,Dn}[, iJt])
    pullback(v::ExteriorVectorValue{K,Dn}, Jt::TensorValue{Dm,Dn}[, iJt])

Pushforward `φ_*ω` of a K-form and pullback `φ^*v` of a K-vector, the two
transports that go against the direction of `φ : ℝᴰᵐ → ℝᴰⁿ`. They are
respectively the pullback and the pushforward by the inverse map `ψ = φ⁻¹`,
whose transposed Jacobian is `iJt = pinvJt(Jt)`, so they are computed from
`iJt` rather than from `Jt` and require `Int(Dm) ≤ Int(Dn)`, see [`pinvJt`](@ref).

Pass a precomputed `iJt` to avoid recomputing it on every call; it is trusted
to equal `pinvJt(Jt)` and is only checked against its expected size.
"""
@generated function pushforward(
  ω::ExteriorFormValue{K,Dw,Tw}, Jt::TensorValue{Dm,Dn,Tj}) where {K,Dw,Dm,Dn,Tw,Tj}

  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  err = _inverse_pullpush_error(Dm, Dn, "pushforward")
  err === nothing || return err
  :( pullback(ω, pinvJt(Jt)) )
end

@generated function pushforward(
  ω::ExteriorFormValue{K,Dw,Tw}, Jt::TensorValue{Dm,Dn,Tj}, iJt::TensorValue{Di1,Di2,Ti}
  ) where {K,Dw,Dm,Dn,Tw,Tj,Di1,Di2,Ti}

  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  err = _pinvJt_shape_error(Dm, Dn, Di1, Di2)
  err === nothing || return err
  err = _inverse_pullpush_error(Dm, Dn, "pushforward")
  err === nothing || return err
  :( pullback(ω, iJt) )
end

@generated function pullback(
  v::ExteriorVectorValue{K,Dw,Tw}, Jt::TensorValue{Dm,Dn,Tj}) where {K,Dw,Dm,Dn,Tw,Tj}

  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  err = _inverse_pullpush_error(Dm, Dn, "pullback")
  err === nothing || return err
  :( pushforward(v, pinvJt(Jt)) )
end

@generated function pullback(
  v::ExteriorVectorValue{K,Dw,Tw}, Jt::TensorValue{Dm,Dn,Tj}, iJt::TensorValue{Di1,Di2,Ti}
  ) where {K,Dw,Dm,Dn,Tw,Tj,Di1,Di2,Ti}

  Dm isa Cov || return :(throw(ArgumentError(
    "the domain (first) index of Jt should be covariant, got size $(size(Jt))"
  )))
  err = _pinvJt_shape_error(Dm, Dn, Di1, Di2)
  err === nothing || return err
  err = _inverse_pullpush_error(Dm, Dn, "pullback")
  err === nothing || return err
  :( pushforward(v, iJt) )
end

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

function flat_inner(a::SkewSymKTensorValue{V,K,D,Ta},
                    b::SkewSymKTensorValue{V,K,D,Tb}) where {V,K,D,Ta,Tb}
  K > D && return zero(Base.promote_op(*,Ta,Tb))
  factorial(K) * flat_alt_inner(a, b)
end

function inner(a::SkewSymKTensorValue{V,K,D,Ta}, b::SkewSymKTensorValue{V,K,D,Tb},
               metric_args::Vararg{_MetricArg}) where {V,K,D,Ta,Tb}
  V === Int && return flat_inner(a, b)
  md = _require_metric_data(a, b, metric_args...)
  K > D && return alt_inner(a, b, md)
  factorial(K) * alt_inner(a, b, md)
end


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


# Christoffel symbols

# Inverse of `_2d_sym_tensor_linear_index` on sorted pairs: the pair `(i,j)`,
# `i ≤ j`. Implemented to be constant foldable.
function _2d_sym_tensor_index_pair(D,l)
  r = l
  for i in 1:D
    r <= D - i + 1 && return (i, i + r - 1)
    r -= D - i + 1
  end
  @unreachable "a $D×$D symmetric tensor stores no component $l"
end

"""
    partials_to_Christoffel_1(∂g, metric_args...)

Christoffel symbols of the first kind `Γ` of a metric `g`, from the
first order partial derivatives `∂g = partial_derivatives(g)(x)` of its value,

    Γₖᵢⱼ = (∂ᵢgⱼₖ + ∂ⱼgᵢₖ - ∂ₖgᵢⱼ) / 2

See also [`partials_to_Christoffel_2`](@ref).
"""
function partials_to_Christoffel_1(
  ∂g::ProductTensorValue{A,SymTensorValue{Dg},T,Tuple{Dk,Dg,Dg},3,L}
) where {A,Dk,Dg,T,L}

  Dg isa Cov || throw(ArgumentError(
    "the Christoffel symbols require a covariant metric tensor, got metric indices of variance $(typeof(Dg))"))
  Dk === Dg || throw(DimensionMismatch(
    "the derivation index and the metric indices must share dimension and variance, got $Dk and $Dg"))
  TT = Base.promote_op(/,T,Int)
  iszero(L) &&
    return zero(ProductTensorValue{CoVectorValue{Dg},SymTensorValue{Dg},TT})

  D = Int(Dg)
  data = ntuple(Val(L)) do n
    @inline
    sij, k = fldmod1(n, D)
    i, j = _2d_sym_tensor_index_pair(D, sij)
    sjk = _2d_sym_tensor_linear_index(D,j,k)
    sik = _2d_sym_tensor_linear_index(D,i,k)
    (indep_comp_getindex(∂g, (sjk-1)*D + i) +
     indep_comp_getindex(∂g, (sik-1)*D + j) -
     indep_comp_getindex(∂g, (sij-1)*D + k)) / 2
  end

  ProductTensorValue{CoVectorValue{Dg},SymTensorValue{Dg}}(data)
end

"""
    partials_to_Christoffel_2(∂g, metric_args...)

Christoffel symbols of the second kind `Γ^k_{i j}` of a metric `g`, from the
first order partial derivatives `∂g = partial_derivatives(g)(x)` of its value,

    Γᵏᵢⱼ = g_invᵏᵐ (∂ᵢgⱼₘ + ∂ⱼgᵢₘ - ∂ₘgᵢⱼ) / 2

`metric_args` contains the metric information, see [`parse_metric_data`](@ref).
`partials_to_Christoffel_2` uses `g_inv`.

See also [`partials_to_Christoffel_1`](@ref)
"""
function partials_to_Christoffel_2(
  ∂g::ProductTensorValue{A,SymTensorValue{Dg},T,Tuple{Dk,Dg,Dg},3,L},
  metric_args::Vararg{_MetricArg}) where {A,Dk,Dg,T,L}

  Γ1 = partials_to_Christoffel_1(∂g)
  V = ProductTensorValue{ContraVectorValue{inverse_variance(Dg)},SymTensorValue{Dg}}
  iszero(L) && return zero(V{eltype(Γ1)})

  md = parse_metric_data(metric_args...)
  isidentity(md) && return V(Tuple(Γ1))

  (; g_inv) = md
  D = Int(Dg)
  _check_metric_dim(g_inv, D)

  # `Γ_{kij} = g_{kl} Γˡᵢⱼ`
  Γ = SMatrix{D,(D*(D+1))÷2}(Tuple(Γ1))
  V(Tuple(get_array(g_inv) * Γ))
end
