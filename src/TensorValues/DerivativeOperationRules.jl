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

for op in (:*,:⋅,:⊙,:⊗,:∧)
  @eval partials_op_rule(::Val{1}, ::typeof($op)) = Base.Fix1(_partials_product_rule, $op)

  @eval partials_op_rule(::Val{2}, ::typeof($op)) =
    Base.Fix1(_second_partials_product_rule, $op)
end
