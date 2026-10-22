# Helpers to implement derivative rules on ProductTensorValue representing the
# partial derivatives of a tensor

# The order `O` partial derivatives, with respect to a `D` dimensional
# coordinate, of a `B` valued tensor.
const _PartialsTensor{O,D,B} = ProductTensorValue{SymNTensorValue{NTuple{O,D}},B}

"""
    _leibniz_cross_term(op, ∂f1, ∂f2)

The mixed term `∂ᵢf₁ op ∂ⱼf₂ + ∂ⱼf₁ op ∂ᵢf₂` of the order 2 Leibniz rule, as a
product tensor holding the two derivation indices in its first factor.

The term is symmetric in `i` and `j` by construction, which is what lets the two
derivation indices share one fully symmetric tensor.
"""
@generated function _leibniz_cross_term(
  op, ∂f1::_PartialsTensor{1,D}, ∂f2::_PartialsTensor{1,D}) where D

  A = SymNTensorValue{Tuple{D,D}}
  comps = [ :( op(as[$i], bs[$j]) + op(as[$j], bs[$i]) )
            for (i,j) in _indep_comp_cartesian_indices(A) ]
  quote
    as, bs = _second_factors(∂f1), _second_factors(∂f2)
    _product_from_second_factors($A, $(Expr(:tuple, comps...)))
  end
end

# Whether a contraction of a product tensor falls inside its second factor
_contracts_in_second_factor(::Type{P}) where P<:ProductTensorValue =
  !iszero(length(size(_second_factor(P)))) && !iszero(num_indep_components(P))

# For `∂_I(ϕ_J) = Σₐ ∂_I(bₐ(x)) ⊗ tₐ_J` where `bₐ` are Fields, and `tₐ` constant MultiValue
function outer(a::_PartialsTensor, b::MultiValue)
  !iszero(num_indep_components(a)) || return contracted_product(Val(0), a, b)
  _map_second_factor(x -> outer(x, b), a)
end

# For `∂_I(ϕ⋅b) = ∂_I(ϕ)⋅b` with `b` constant
function dot(a::_PartialsTensor, b::MultiValue)
  _contracts_in_second_factor(typeof(a)) || return contracted_product(Val(1), a, b)
  _map_second_factor(x -> dot(x, b), a)
end

@generated function dot(a::ProductTensorValue{Aa,Ba}, b::ProductTensorValue{Ab,Bb}) where {Aa,Ba,Ab,Bb}

  Sa, Sb = size(Ba), size(Ab)
  if length(Sa) != 1 || length(Sb) != 1
    return :( @notimplemented "A product tensor contraction consumes a whole factor
      of each operand, so it needs a single index in the second factor of the first
      operand and in the first factor of the second one, got $($Ba) and $($Ab)" )
  end

  i, j = only(Sa), only(Sb)
  can_contract(i, j) || return :( throw(DimensionMismatch(
    "Indices $($i) and $($j) cannot be contracted, see [`can_contract`](@ref)")) )

  iLa, Dc = num_indep_components(Aa), Int(i)
  comps = [ :( +($([:( as[$la][$l] * bs[$l] ) for l in 1:Dc]...)) ) for la in 1:iLa ]
  quote
    as, bs = _second_factors(a), _second_factors(b)
    _product_from_second_factors($Aa, $(Expr(:tuple, comps...)))
  end
end

# fixes method ambiguity of the two above
function dot(a::_PartialsTensor, b::ProductTensorValue)
  # Whether a product tensor contraction consumes a whole factor of each operand,
  # which it does when a single index is left to contract on either side.
  contracts_whole_factors = isone(length(size(_second_factor(a)))) && isone(length(size(_first_factor(b))))

  contracts_whole_factors  ?
    @invoke(dot(a::ProductTensorValue, b::ProductTensorValue)) :
    @invoke(dot(a::_PartialsTensor,    b::MultiValue))
end


#######################################
# API for differential operator rules #
#######################################

"""
    partials_op_rule(Val(N), op) -> op_rule
    partials_op_rule(op) -> op_rule

The differentiation rule used to commute an operation `op` with the `N`th order
partial derivative operator, that is

    x -> PartialDerivatives{N}()(op(args(x)...))

can be re-written as

    x -> op_rule(args(x)...[, ∂args(x)...][, ...], ∂ᴺargs(x)...) # if given N
    x -> op_rule(∂ᴺargs(x)...) # if no N

where `∂args` are the `partial_derivatives` of the `args`, `∂∂args` are
their `second_partial_derivatives`, etc. When `Val(N)` is given, `op_rule`
assumes that all partial derivatives of all arguments until derivative order
`N` are given.

If `Val(N)` is omitted, the rule should be valid for any `N`, which requires
`op` to be linear jointly in all its argument. In this case, the argument is
assumed to be the `N`th derivatives of the arguments only.
"""
function partials_op_rule end

# jointly linear operations

partials_op_rule(::typeof(+)) = +
partials_op_rule(::typeof(-)) = -
partials_op_rule(::typeof(flat_hodge_star)) = Base.Fix1(_map_second_factor, flat_hodge_star)


# bilinear product rules, order 1 and 2 only

"""
    _partials_product_rule(op, f1, f2, ∂f1, ∂f2)

Leibniz rule for a bilinear `op`: it acts on the second factor of each argument's
partial derivatives against the other argument's value, the first factor holding
the derivation index.
"""
_partials_product_rule(op, f1, f2, ∂f1, ∂f2) =
  _map_second_factor(b -> op(b,f2), ∂f1) + _map_second_factor(b -> op(f1,b), ∂f2)

"""
    _second_partials_product_rule(op, f1, f2, ∂f1, ∂f2, ∂∂f1, ∂∂f2)

Order 2 Leibniz rule for a bilinear `op`: each argument's second derivatives act
against the other argument's value, and the mixed term pairs their first
derivatives.
"""
_second_partials_product_rule(op, f1, f2, ∂f1, ∂f2, ∂∂f1, ∂∂f2) =
  _map_second_factor(b -> op(b,f2), ∂∂f1) +
  _leibniz_cross_term(op, ∂f1, ∂f2) +
  _map_second_factor(b -> op(f1,b), ∂∂f2)

for op in (:*,:⋅,:⊙,:⊙ᴱ,:⊗,:∧,:ι,:koszul,:alt_endomorphism)
  @eval partials_op_rule(::Val{1}, ::typeof($op)) = Base.Fix1(_partials_product_rule, $op)

  @eval partials_op_rule(::Val{2}, ::typeof($op)) =
    Base.Fix1(_second_partials_product_rule, $op)
end


# metric dependent operations, which are not linear in the metric
# TODO review that and look for simplification

"""
    _partials_inv_rule(a, ∂a)

Order 1 rule for `inv`: `∂_c(a⁻¹) = -a⁻¹ (∂_c a) a⁻¹`.
"""
@inline _partials_inv_rule(a::Number, ∂a) = _map_second_factor(b -> -b/a^2, ∂a)

@inline function _partials_inv_rule(a::MultiValue, ∂a)
  a_inv = inv(a)
  _map_second_factor(b -> -a_inv⋅b⋅a_inv, ∂a)
end

# `a⁻¹ (∂_c a) a⁻¹` is a congruence of a symmetric tensor, hence symmetric
@inline function _partials_inv_rule(a::SymTensorValue, ∂a)
  a_inv = inv(a)
  _map_second_factor(b -> -congruent_prod(b, a_inv), ∂a)
end

"""
    _partials_sqrt_abs_det_rule(g, ∂g)

Order 1 rule for `_sqrt_abs_det`, Jacobi's formula:

    ∂_c √|det g| = √|det g| g^{ab} (∂_c g_{ab}) / 2
"""
@inline function _partials_sqrt_abs_det_rule(g::SymTensorValue, ∂g)
  half_sqrt_abs_g = _sqrt_abs_det(g) / 2
  g_inv = inv(g)
  _map_second_factor(b -> half_sqrt_abs_g * double_contraction(g_inv, b), ∂g)
end

# Σ_r of the determinant of the `K`×`K` submatrix of `M` with rows `I` and
# columns `J`, with its row `r` taken from `∂M`: the derivative of that minor,
# which is `K`-linear in its rows.
function _minor_derivative_expr(M::Symbol, ∂M::Symbol, I, J)
  K = length(I)
  iszero(K) && return :(zero(eltype($M)))
  Expr(:call, :+, (_det_expr((a,b) -> :($(a == r ? ∂M : M)[$(I[a]), $(J[b])]), K)
                   for r in 1:K)...)
end

# `_alt_compound_comps` differentiated, the `K`-th compound matrix of `M` being
# `K`-linear in the rows of its minors.
_alt_compound_derivative_comps(M::Symbol, ∂M::Symbol, c_out, c_in, T, src::Symbol) =
  [ _sum_expr([ :( $(_minor_derivative_expr(M, ∂M, I, J)) * indep_comp_getindex($src, $l) )
                for (l, J) in enumerate(c_in) ], T)
    for I in c_out ]

#     _hodge_star_slice(ω, ∂ωc, M, ∂Mc, vol, dvolc)
#
# The derivation slice `∂_c(⋆ω)` of the derivative of [`hodge_star`](@ref), from
# the `c`-th slices `∂ωc`, `∂Mc` and `dvolc` of the derivatives of its arguments.
#
# `⋆ω = vol ⋆ᴱ(Mᴷω)` raises the `K` indices of `ω` through the `K`-th compound of
# the metric `M` of variance opposite to `ω`, so the rule is a product rule over
# `vol` and the `K` rows of each minor of that compound.
@generated function _hodge_star_slice(
  ω::SkewSymKTensorValue{V,K,Dw,Tw}, ∂ωc::SkewSymKTensorValue{V,K,Dw},
  M::SymTensorValue, ∂Mc::MultiValue, vol::Number, dvolc::Number) where {V,K,Dw,Tw}

  D = Int(Dw)
  K <= D || return :(@unreachable $("hodge star requires K ≤ D, got K = $K and D = $D"))
  K <= 4 || return :(@notimplemented $("the derivative of a Hodge star expands the \
  Leibniz sum of a K×K determinant K times over, which is $(factorial(K)*K) terms \
  per minor at K = $K"))

  Kc  = D - K
  L   = binomial(D, K)
  c_K = sorted_combinations(D, K)
  T   = Base.promote_op(*, Tw, eltype(M))

  swapped    = _alt_compound_comps(:M, c_K, c_K, T, :ω)
  d_swapped  = _alt_compound_comps(:M, c_K, c_K, T, :∂ωc)
  dM_swapped = _alt_compound_derivative_comps(:M, :∂Mc, c_K, c_K, T, :ω)

  s = _hodge_star_signs(K, D)
  comps = [_signed(s[L+1-m], :( dvolc * ω_swapped[$(L+1-m)] +
                                vol * (∂ω_swapped[$(L+1-m)] + ∂M_swapped[$(L+1-m)]) ))
           for m in 1:L]

  quote
    ω_swapped  = ($(swapped...),)
    ∂ω_swapped = ($(d_swapped...),)
    ∂M_swapped = ($(dM_swapped...),)
    SkewSymKTensorValue{$V,$Kc,$Dw}(($(comps...),))
  end
end

#     _partials_hodge_star_rule(ω, metric_args..., ∂ω, ∂metric_args...)
#
# Order 1 rule for [`hodge_star`](@ref), one [`_hodge_star_slice`](@ref) per
# derivation index.
#
# Both the metric tensor of the variance `⋆` raises the indices of `ω` with and
# the volume factor `√|det g|` must be among `metric_args`: the rule
# differentiates the arguments it is handed, and derives neither of the two from
# the other. `get_metric_data` returns all of them.
@generated function _partials_hodge_star_rule(args::Vararg{Any,N}) where N

  nm = N ÷ 2 - 1
  Tω = args[1]
  Tω <: SkewSymKTensorValue || return :(@notimplemented $("`hodge_star` is \
  differentiated on an alternating tensor with variance, got $Tω"))
  V, Kw, Dw = Tω.parameters[1], Tω.parameters[2], Tω.parameters[3]

  # The metric of opposite variance to ω is the one `hodge_star` reads, and the
  # volume factor is the sole scalar amongst the metric arguments.
  Dm = inverse_variance(Dw)
  im = findfirst(i -> args[1+i] <: SymTensorValue{Dm}, 1:nm)
  iv = findfirst(i -> args[1+i] <: Number && !(args[1+i] <: MultiValue), 1:nm)
  (isnothing(im) || isnothing(iv)) && return :(@notimplemented $("the derivative \
  of `hodge_star` of a $Kw-tensor of variance $V reads the metric tensor of \
  variance $(typeof(Dm)) and the volume factor √|det g| off its arguments, and \
  derives neither of the two from the other: pass both of them"))

  D = Int(Dw)
  A = _first_factor(args[2+nm])
  pm, pv = 1+im, 1+iv            # the metric tensor and the volume factor
  qm, qv = 2+nm+im, 2+nm+iv      # their derivatives

  # Raising the indices of a form pairs with √|det(g)|, lowering those of a
  # K-vector with its reciprocal, as in `hodge_star` itself.
  _vol  = V === Cov ? :( args[$pv] ) : :( inv(args[$pv]) )
  _dvol(c) = V === Cov ? :( indep_comp_getindex(args[$qv], $c) ) :
                         :( -indep_comp_getindex(args[$qv], $c) / args[$pv]^2 )
  slices = [:( _hodge_star_slice(
                 args[1], ∂ω_slices[$c], M, ∂M_slices[$c], vol, $(_dvol(c))) )
            for c in 1:D]

  quote
    M = args[$pm]
    _check_metric_dim(M, $D)
    ∂ω_slices, ∂M_slices = _second_factors(args[$(2+nm)]), _second_factors(args[$qm])
    vol = $_vol
    _product_from_second_factors($A, ($(slices...),))
  end
end

partials_op_rule(::Val{1}, ::typeof(inv))            = _partials_inv_rule
partials_op_rule(::Val{1}, ::typeof(_sqrt_abs_det))  = _partials_sqrt_abs_det_rule
partials_op_rule(::Val{1}, ::typeof(hodge_star))     = _partials_hodge_star_rule
