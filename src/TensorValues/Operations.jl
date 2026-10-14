###############################################################
# Comparison
###############################################################

(==)(a::Number, b::MultiValue) = false
(==)(a::MultiValue, b::MultiValue) = false
(==)(a::MultiValue{S}, b::MultiValue{S}) where S = a.data == b.data

function (≈)(a::MultiValue,b::MultiValue;kwargs...)
  Sa, Sb = size(a), size(b)
  Sa !== Sb && throw(DimensionMismatch("Tensors of size $(size(a)) and $(size(b)) are not comparable."))
  ≈(get_array(a),get_array(b);kwargs...)
end

function (≈)(
  a::AbstractArray{<:MultiValue}, b::AbstractArray{<:MultiValue}; kwargs...)
  if size(a) != size(b); return false; end
  for (ai,bi) in zip(a,b)
    if !≈(ai,bi;kwargs...); return false; end
  end
  true
end

function isless(a::MultiValue{Tuple{L}}, b::MultiValue{Tuple{L}}) where {L}
  for d in L:-1:1
    if isless(a[d], b[d])
      return true
    elseif isless(b[d], a[d])
      return false
    else
      continue
    end
  end
  false
end

# Comparison uses get_array for robustness to tensors with symmetry or tracelessness
isless(a::Number, b::MultiValue) = all(isless(a, bi) for bi in get_array(b))
function <=(a::Number, b::MultiValue)
  all(a <= bi for bi in get_array(b)) # default doesen't work because a==b is always false
end
isless(a::MultiValue, b::MultiValue) = @unreachable "Comparison is not defined between tensor of order greater than 1"
function <=(a::MultiValue, b::MultiValue)
  isless(a,b) || isequal(a,b)
end

###############################################################
# promotion and conversions
###############################################################

# No default promotion /conversion for tensors of different type names
promote_rule(::Type{<:MultiValue}, ::Type{<:MultiValue}) = Union{}

# But promotion and conversion between the different types of square tensors.
promote_rule(::Type{<:TensorValue{D,D,Ta}},     ::Type{<:SymTensorValue{D,Tb}})          where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:TensorValue{D,D,Ta}},     ::Type{<:SkewSymTensorValue{D,Tb}})      where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:TensorValue{D,D,Ta}},     ::Type{<:SymTracelessTensorValue{D,Tb}}) where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:SymTensorValue{D,Ta}},    ::Type{<:SkewSymTensorValue{D,Tb}})      where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:SymTensorValue{D,Ta}},    ::Type{<:SymTracelessTensorValue{D,Tb}}) where {D,Ta,Tb} = SymTensorValue{D,promote_type(Ta,Tb)}
promote_rule(::Type{<:SkewSymTensorValue{D,Ta}},::Type{<:SymTracelessTensorValue{D,Tb}}) where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}

convert(::Type{<:TensorValue{D,D,Ta}},  a::SymTensorValue{D,Tb})          where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}(get_array(a))
convert(::Type{<:TensorValue{D,D,Ta}},  a::SkewSymTensorValue{D,Tb})      where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}(get_array(a))
convert(::Type{<:TensorValue{D,D,Ta}},  a::SymTracelessTensorValue{D,Tb}) where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}(get_array(a))
convert(::Type{<:SymTensorValue{D,Ta}}, a::SymTracelessTensorValue{D,Tb}) where {D,Ta,Tb} = SymTensorValue{D,promote_type(Ta,Tb)}(a.data)

# TODO deprecate next two methods ? A few stuff depend on this behavior but very ugly.
function convert(V::Type{<:MultiValue}, a::T) where T<:_Scalar
  isone(length(V)) && isone(num_indep_components(V)) || error("Cannot convert value of type $V to type $T")
  V(a)
end
function convert(T::Type{<:_Scalar}, a::V) where V<:MultiValue
  isone(length(a)) || error("Cannot convert value of type $V to type $T")
  T(a[1])
end

###############################################################
# Addition / subtraction
###############################################################

Base.iszero(a::MultiValue) = all(iszero.(a.data))

for op in (:+,:-)
  @eval begin

    function ($op)(a::T) where T<:MultiValue
      r = map($op, get_indep_components(a))
      T(r)
    end

    function ($op)(a::V, b::V) where V<:MultiValue
      r = map(($op), get_indep_components(a), get_indep_components(b))
      V(r)
    end
  end
end


###############################################################
# Arithmetics with UniformScaling
###############################################################

# Dealing with LinearAlgebra.I fallbacks due to `MultiValue <: Number`

const _UniformScalingMsg = "UniformScaling can only be added to or subtracted \
  from a square order-2 tensor, got size "

for op in (:+,:-)
  @eval begin
    function ($op)(a::MultiValue{Tuple{D1,D2}}, ::UniformScaling) where {D1,D2}
      Int(D1) == Int(D2) || throw(DimensionMismatch(_UniformScalingMsg*"$(size(a))."))
      ($op)(a, one(a))
    end

    function ($op)(::UniformScaling, a::MultiValue{Tuple{D1,D2}}) where {D1,D2}
      Int(D1) == Int(D2) || throw(DimensionMismatch(_UniformScalingMsg*"$(size(a))."))
      ($op)(one(a), a)
    end

    function ($op)(a::MultiValue, ::UniformScaling)
      throw(DimensionMismatch(_UniformScalingMsg*"$(size(a))."))
    end

    function ($op)(::UniformScaling, a::MultiValue)
      throw(DimensionMismatch(_UniformScalingMsg*"$(size(a))."))
    end
  end
end

(*)(a::MultiValue, J::UniformScaling) = a * J.λ
(*)(J::UniformScaling, a::MultiValue) = J.λ * a
(/)(a::MultiValue, J::UniformScaling) =
  iszero(J.λ) ? throw(SingularException(1)) : a / J.λ
(\)(J::UniformScaling, a::MultiValue) =
  iszero(J.λ) ? throw(SingularException(1)) : a / J.λ

###############################################################
# Matrix Division
###############################################################

# Sum of `terms`, without the wrapping call when there is only one of them.
_plus_expr(terms) = length(terms) == 1 ? only(terms) : Expr(:call, :+, terms...)

# `Σₗ wₗ * term(l)` over the stored components `ls` of a tensor layout weighted
# by their multiplicities `ws`, see [`_stored_comp_weights`](@ref): the stored
# components sharing a weight are summed first, so each weight
# multiplies once. `ls` is not empty.
function _weighted_sum_expr(term, ls, ws)
  groups = Dict{Int,Vector{Any}}()
  for (l, w) in zip(ls, ws)
    push!(get!(Vector{Any}, groups, w), term(l))
  end
  _plus_expr([ isone(w) ? _plus_expr(ts) : :( $w * $(_plus_expr(ts)) )
               for (w, ts) in sort!(collect(groups), by=first) ])
end

"""
    (\\)(a::MultiValue{Tuple{D1,D2}}, b::MultiValue)

Left division `inv(a)⋅b`, that is the solution `x` of `a⋅x = b`.

`a` must be square, `D1` and `D2` being equal in value. The contracted pair is
the first index of `inv(a)` -- the second index of `a` with flipped variance --
against the first index of `b`, so the first indices of `a` and `b` must be
identical, variance included. The result carries the indices
`(inverse_variance(D2), size(b)[2:end]...)`.
"""
@generated function (\)(a::MultiValue{Tuple{D1,D2}}, b::MultiValue{Sb,Tb,Nb}) where {D1,D2,Sb,Tb,Nb}
  if iszero(Nb)
    msg = "Left division by an order-0 tensor is not implemented, got size $Sb"
    return :( @notimplemented $msg )
  end
  if Int(D1) != Int(D2)
    msg = "Left division requires a square left operand, got dimensions $D1 and $D2"
    return :( throw(DimensionMismatch($msg)) )
  end
  B1 = first(size(b))
  if D1 !== B1
    msg = "Cannot left-divide a tensor of first index $B1 by a tensor of first
      index $D1, they must be identical, variance included."
    return :( throw(DimensionMismatch($msg)) )
  end

  # `inv(a)⋅b` is not symmetric even when `b` is, so the result type is built
  # from the indices and never from `typeof(b)`.
  Sr = (inverse_variance(D2), size(b)[2:Nb]...)
  V = _tensor_base_type_expr(Nb, Sr)
  quote
    r = get_array(a) \ get_array(b)
    $V{eltype(r)}(Tuple(r))
  end
end

"""
    (/)(a::MultiValue{Tuple{Da1,Da2}}, b::MultiValue{Tuple{Db1,Db2}})

Right division `a⋅inv(b)`, that is the solution `x` of `x⋅b = a`.

`b` must be square, `Db1` and `Db2` being equal in value. The contracted pair is
the second index of `a` with the second index of `inv(b)` -- which is `Db1` with
flipped variance -- so the second indices of `a` and `b` must be identical,
variance included. The result carries the indices `(Da1, inverse_variance(Db1))`.

Those two conditions are the ones `\\` enforces on the adjoints, which is how
this is assembled, following `LinearAlgebra`'s generic `A/B = (Bᵀ \\ Aᵀ)ᵀ`.
"""
(/)(a::MultiValue{Tuple{Da1,Da2}}, b::MultiValue{Tuple{Db1,Db2}}) where {Da1,Da2,Db1,Db2} =
  adjoint(adjoint(b) \ adjoint(a))

###############################################################
# Operations with scalar numbers
###############################################################

@generated function _bc(f, a::NTuple{N}, b::_Scalar) where N
  res = Expr(:tuple)
  for i in 1:N
    push!(res.args, :(f(a[$i],b)) )
  end
  res
end

@generated function _bc(f, b::_Scalar, a::NTuple{N}) where N
  res = Expr(:tuple)
  for i in 1:N
    push!(res.args, :(f(b,a[$i])) )
  end
  res
end

for op in (:+,:-,:*)
  @eval begin
    function ($op)(a::MultiValue, b::_Scalar)
      r = _bc($op, get_indep_components(a), b)
      T = _eltype($op, r, a, b)
      M = change_eltype(a, T)
      M(r)
    end

    function ($op)(a::_Scalar, b::MultiValue)
      r = _bc($op, a, get_indep_components(b))
      T = _eltype($op, r, a, b)
      M = change_eltype(b, T)
      M(r)
    end
  end
end

const _AbstractTracelessTensor{D} = Union{SymTracelessTensorValue{D},SkewSymTensorValue{D}}
_err = "This operation is undefined for traceless tensors"
(+)(::_AbstractTracelessTensor, ::_Scalar) = error(_err)
(+)(::_Scalar, ::_AbstractTracelessTensor) = error(_err)
(-)(::_AbstractTracelessTensor, ::_Scalar) = error(_err)
(-)(::_Scalar, ::_AbstractTracelessTensor) = error(_err)

function (/)(a::MultiValue, b::_Scalar)
  r = _bc(/, get_indep_components(a), b)
  T = _eltype(/, r, a, b)
  P = change_eltype(a, T)
  P(r)
end

@inline function _eltype(op, r, a)
  eltype(r)
end

@inline function _eltype(op, r, a, b)
  eltype(r)
end

@inline function _eltype(op, r, a...)
  eltype(r)
end

@inline function _eltype(op, r::Tuple{}, a)
  typeof(op(zero(eltype(a))))
end

@inline function _eltype(op, r::Tuple{}, a, b)
  typeof(op(zero(eltype(a)), zero(eltype(b))))
end

@inline function _eltype(op, r::Tuple{}, a...)
  typeof(reduce(op, zero.(eltype.(a))))
end


###############################################################
# Scalar product
###############################################################

function (*)(a::MultiValue, b::MultiValue)
  msg = """
  Method (*)(::$(typeof(a)),::$(typeof(b))) has been removed.
  Depending the case, use simple contraction dot aka ⋅ (\\cdot) or full contraction inner aka ⊙ (\\odot) instead.
  """
  error(msg)
end

# Authorise * for order 0 tensors, they have no index to contract but represent a scalar component.
(*)(a::MultiValue{Tuple{}}, b::MultiValue) = a[1]*b
(*)(a::MultiValue, b::MultiValue{Tuple{}}) = a*b[1]
(*)(a::MultiValue{Tuple{}}, b::MultiValue{Tuple{}}) = ScalarValue(a[1]*b[1])
# Of two order 0 factors, a `ScalarValue` holds no shape of its own, so the other
# one's type is kept; two forms of one family keep their embedding dimension.
(*)(a::ScalarValue, b::MultiValue{Tuple{}}) = a[1]*b
(*)(a::MultiValue{Tuple{}}, b::ScalarValue) = a*b[1]
(*)(a::ScalarValue, b::ScalarValue) = ScalarValue(a[1]*b[1])
(*)(a::SkewSymKTensorValue{V,0,D}, b::SkewSymKTensorValue{V,0,D}) where {V,D} = a*b[1]

###############################################################
# Dot product (simple contraction)
###############################################################

"""
    dot(a::MultiValue{Tuple{...,D}}, b::MultiValue{Tuple{D,...}})
    a ⋅¹ b
    a ⋅ b

Single contraction of two tensors `a` and `b`, of the last index of `a` with
the first index of `b`. The corresponding dimensions `D` must match. No symmetry
is preserved.
On two vectors, this is the same as the inner product.

A `VectorValue`, whose index is un-varianced, is assumed contravariant: it
also contracts against a covariant index, and the indices left over keep
their own variance, so a vector contracted with a covariant two-tensor gives
a covector.
"""
dot(a::MultiValue, b::MultiValue) = contracted_product(Val(1), a, b)

const ⋅¹ = dot

###############################################################
# Inner product (full contraction)
###############################################################

inner(a::_Scalar, b::_Scalar, metric_args::Vararg{_MetricArg}) = a * b

# `Σ_I a_I b_I` read off every component of the two full tensors.
function _dense_flat_inner(
  a::MultiValue{Sa,Ta,N}, b::MultiValue{Sb,Tb,N}) where {Sa,Ta,Sb,Tb,N}
  iszero(length(a)) && return zero(Base.promote_op(*,Ta,Tb))
  # This is StaticArrays's dot, but that doesn't conjugate a.
  za, zb = zero(Ta), zero(Tb)
  ret = za*zb + za*zb
  @inbounds @simd for j in 1:length(a)
    ret += a[j]*b[j]
  end
  ret
end

# The body of the Euclidean contraction of two tensors `a` and `b`,
#
# For tensors with same stored components, so the sum over the full tensors is
# the sum over the stored components weighted by their multiplicities.
@generated function _flat_inner_body(a, b)
  dense_e = :( _dense_flat_inner(a, b) )
  a_sci, b_sci = try
    stored_comp_indices(a), stored_comp_indices(b)
  catch
    return dense_e
  end
  get_array(a_sci) == get_array(b_sci) || return dense_e
  num_stored_components(a) < length(a) || return dense_e

  ws = _stored_comp_weights(a) # shouln't error if stored_comp_indices didn't
  isempty(ws) && return dense_e
  quote
    da, db = Tuple(a), Tuple(b)
    $(_weighted_sum_expr(k -> :( da[$k]*db[$k] ), eachindex(ws), ws))
  end
end

"""
    flat_inner(ω::MultiValue{S}, η::MultiValue{S})
    ⊙ᴱ(ω,η)
    ω ⊙ᴱ η

Tensor inner product in Euclidean space, that is `inner`/`⊙` even if `ω` and
`η` have variance but with iddentity (inverse) metric:

    ⊙ᴱ(ω,η) = Σ_I ω_I η_I

over every cartesian index `I`.
"""
function flat_inner(a::MultiValue{S,Ta,N}, b::MultiValue{S,Tb,N}) where {S,Ta,Tb,N}
  iszero(length(a)) && return zero(Base.promote_op(*,Ta,Tb))
  _flat_inner_body(a, b)
end

const ⊙ᴱ = flat_inner

"""
    inner(a::MultiValue{S}, b::MultiValue{S}, metric_args...) -> scalar
    inner(a, b, LinearAlgebra.I) = flat_inner(a,b) = a ⊙ᴱ b
    a ⊙ b

Inner product of two tensors, that is the inner product along each indices. `a`
and `b` must have iddentical sizes, including variances: the two live in the
same tensor space, which is the only case in which pairing index `k` of `a`
with index `k` of `b` is an invariant operation.

The metric information is strictly required for tensors with variance since the
inner product requires lowering and or raising the varianced indices of one
argument, see the [metric operation](@ref "Metric dependent operations")
documentation. `metric_args` are the arguments of
[`parse_metric_data`](@ref), in any order.

The indices without variance are contracted assuming iddentity metric. If `a`
and `b` are type `(r,s,u)` tensors, the mathematical definition of inner is
thus

    a ⊙ b = Σ_{IJKLM} a^{i1, ..., ir}_{j1, ..., js}{k1, ..., ku} *
              g_{i1,l1} * ... * g_{ir,lr} *
              ginv^{j1,m1} * ... * ginv^{js,ms} *
            b^{l1, ..., lr}_{m1, ..., ms}{k1, ..., ku}

where the sum is over all the cartesian indices `I` and `L` of the `r`
contravariant indices, `J` and `M` that on the `s` covariant indices, and `K`
of the `u` un-varianced indices. A contravariant index pair is thus contracted
through `g` and a covariant one through `g_inv`, as in [`alt_inner`](@ref).

The Euclidean method [`flat_inner`](@ref) is called if the given metric is
`LinearAlgebra.I`.

# Extended help

Neither argument is conjugated: `inner` is bilinear, so for complex components
`a ⊙ b = Σ_I a_I b_I` and `a ⊙ a` is not `|a|²`. The Galerkin weak forms of
time-harmonic problems are complex-symmetric rather than Hermitian, so they are
written directly as `a ⊙ b`; conjugate an argument explicitly, as
`inner(a, conj(b))`, for the Hermitian product -- which is what [`norm`](@ref)
does internally.
"""
function inner(a::MultiValue, b::MultiValue, metric_args::Vararg{_MetricArg})
  @notimplemented "Sizes of tensors must match."
end

function inner(a::MultiValue{S,Ta,N}, b::MultiValue{S,Tb,N},
                          metric_args::Vararg{_MetricArg}) where {S,Ta,Tb,N}

  sa = size(a)

  # Call flat_inner if S has no variance, md is ignored
  all(i -> i isa Int, sa) && return flat_inner(a, b)
  md = _require_metric_data(a, b, metric_args...)

  # A tensor with no component contracts to zero whatever the metric.
  iszero(length(a)) && return zero(Base.promote_op(*,Ta,Tb))
  isidentity(md) && return flat_inner(a, b)

  # `flat` and `sharp` act on all the indices they are given at once, and check
  # both their variance and the dimension of the metric against them.
  indices = ntuple(identity, Val(N))
  cons = filter(k -> sa[k] isa Con, indices)
  covs = filter(k -> sa[k] isa Cov, indices)
  b♯  = isempty(cons) ? b  : flat( b, Val(cons), md; drop_sym=Val(true))
  b♭♯ = isempty(covs) ? b♯ : sharp(b♯,Val(covs), md; drop_sym=Val(true))

  @check all(k -> can_contract(sa[k], size(b♭♯)[k]), 1:N)
  iszero(length(a)) && return zero(Base.promote_op(*,Ta,Tb))
  # this is a full contracted_product optimized for dependent component inputs
  _flat_inner_body(a, b♭♯)
end

# skew symmetric and symmetric tensors are orthogonal
function inner(::SkewSymTensorValue{D,Ta}, ::AbstractSymTensorValue{D,Tb}, ::Vararg{_MetricArg}) where {D,Ta,Tb}
  zero(Base.promote_op(*,Ta,Tb))
end
function inner(::AbstractSymTensorValue{D,Tb}, ::SkewSymTensorValue{D,Ta}, ::Vararg{_MetricArg}) where {D,Ta,Tb}
  zero(Base.promote_op(*,Ta,Tb))
end
function flat_inner(::SkewSymTensorValue{D,Ta}, ::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  zero(Base.promote_op(*,Ta,Tb))
end
function flat_inner(::AbstractSymTensorValue{D,Ta}, ::SkewSymTensorValue{D,Tb}) where {D,Ta,Tb}
  zero(Base.promote_op(*,Ta,Tb))
end

# inner should only be defined on tensors with same size
@deprecate inner(a::MultiValue{Tuple{D,D,D,D}}, b::MultiValue{Tuple{D,D}}) where D double_contraction(a,b)
@deprecate inner(a::MultiValue{Tuple{D,D}}, b::MultiValue{Tuple{D,D,D,D}}) where D double_contraction(a,b)

# The Euclidean multiplicities `wᵢ = eᵢ⊙ᴱeᵢ` of the component basis `{eᵢ}` of
# `V`, or `nothing` if that basis is not orthogonal -- in which case the
# contraction is not a weighted sum over the independent components.
function _component_basis_weights(V::Type{<:MultiValue})
  basis = component_basis(V)
  for i in eachindex(basis), j in eachindex(basis)
    i == j && continue
    iszero(_dense_flat_inner(basis[i],basis[j])) || return nothing
  end
  Int[ _dense_flat_inner(b,b) for b in basis ]
end

# Σ_(la,lb) wa(la)*wb(lb)*a_(la,lb)*b_(la,lb), the component basis of a product
# being the products of the factors' component bases.
@generated function flat_inner(
  a::ProductTensorValue{A,B,Ta,S,N,L}, b::ProductTensorValue{A,B,Tb,S,N,L}) where {A,B,Ta,Tb,S,N,L}

  zero_result = :( zero(Base.promote_op(*,Ta,Tb)) )
  iszero(L) && return zero_result

  generic_e = :( @invoke flat_inner(a::MultiValue{S,Ta,N}, b::MultiValue{S,Tb,N}) )

  # The multiplicities of a component basis are Euclidean, so they are the same
  # whether or not the factors carry variance.
  wa = _component_basis_weights(make_concretetype(A{Int}))
  wb = _component_basis_weights(make_concretetype(B{Int}))
  # A factor whose component basis is not orthogonal has no such weighted sum,
  # so the contraction is read off every component instead.
  (isnothing(wa) || isnothing(wb)) && return generic_e

  La, = _prod_counts(A,B)
  terms = Expr[]
  for lb in eachindex(wb), la in eachindex(wa)
    w = wa[la]*wb[lb]
    iszero(w) && continue
    k = (lb-1)*La + la
    prod_expr = :( a.data[$k]*b.data[$k] )
    push!(terms, isone(w) ? prod_expr : :( $w * $prod_expr ))
  end
  isempty(terms) && return zero_result
  Expr(:call, :+, terms...)
end

const ⊙ = inner

###############################################################
# Double Contractions w/ products
###############################################################

"""
    double_contraction(a::MultiValue{Tuple{...,Da,Ea}}, b::MultiValue{Tuple{Db,Eb,...})
    a ⋅² b

Double contraction of two tensors `a` and `b`, along the two last indices of `a`
and two first of `b`, that is `Da` contracts with `Db` and `Ea` with `Eb`.

The `double_contraction` between second- and/or fourth-order symmetric tensors
preserves the symmetry (returns a symmetric tensor type).
"""
double_contraction(a::MultiValue, b::MultiValue) = contracted_product(Val(2), a, b)

function _comp_prod_double_symfourth4(::Val{D},a,b,i,j,k,l) where D
  s = zero(promote_type(eltype(a), eltype(b)))
  @inbounds for m in 1:D
    for n in 1:D
      s += a[i,j,m,n]*b[m,n,k,l]
    end
  end
  s
end

# c_ijkl = a_ijmn*b_mnkl (general case)
@generated function double_contraction(a::SymFourthOrderTensorValue{D,Ta}, b::SymFourthOrderTensorValue{D,Tb}) where {D,Ta,Tb}
  D isa Int || return :( contracted_product(Val(2), a, b) )
  iszero(D) && return :(SymFourthOrderTensorValue{0,$(Base.promote_op(*,Ta,Tb))}())

  comps = Expr[]
  for j in 1:D
    for i in j:D
      for l in 1:D
        for k in l:D
          if D < 4
            comp = Expr(:call, :+)
            for m in 1:D
              for n in 1:D
                push!(comp.args, :(a[$i,$j,$m,$n]*b[$m,$n,$k,$l]) )
              end
            end
            push!(comps, comp)
          else
            # for D=4, this compiles in <0.05 s, while the other takes 200 s (Julia 1.12).
            # But runtime is 5 times slower 1.6 μs, vs 300 ns
            #
            # I tried passing i,j,k,l by Val() instead of value, compilation
            # takes 1s, the runtime isn't improved
            #
            # This means that the acceleration we get by removing the function
            # likely comes from llvm optimization that recycle partial results
            # of the products for different i,j,k,l
            push!(comps, :(_comp_prod_double_symfourth4(Val(D),a,b,$i,$j,$k,$l)) )
          end
        end
      end
    end
  end
  :( SymFourthOrderTensorValue{D}($(comps...)) )
end

# c_ijk = a_ilm*b_lmjk
@generated function double_contraction(a::ThirdOrderTensorValue{D1,Da,Da,Ta}, b::SymFourthOrderTensorValue{Db,Tb}) where {D1,Da,Db,Ta,Tb}
  can_contract(Da,Db) || return :( contracted_product(Val(2), a, b) ) # for error msg
  D = Int(Da)
  iszero(length(a)) && return :(zero(ThirdOrderTensorValue{D1,Db,Db,$(Base.promote_op(*,Ta,Tb))}))
  comps = Expr[]
  for k in 1:D
    for j in 1:D
      for i in 1:Int(D1)
        summands = [:(a[$i,$l,$m]*b[$l,$m,$j,$k]) for l in 1:D for m in 1:D]
        push!(comps, Expr(:call, :+, summands...))
      end
    end
  end
  :( ThirdOrderTensorValue{D1,Db,Db}($(comps...)) )
end

# c_ij = a_ijkl*b_kl
@generated function double_contraction(a::SymFourthOrderTensorValue{Da,Ta}, b::AbstractSymTensorValue{Db,Tb}) where {Da,Db,Ta,Tb}
  can_contract(Da,Db) || return :( contracted_product(Val(2), a, b) ) # for error msg
  D = Int(Da)
  iszero(D) && return :(zero(SymTensorValue{Da,$(Base.promote_op(*,Ta,Tb))}))
  # `b` is contracted over the whole tensor, read off the stored comp indices
  # of the symmetric layout weighted by their multiplicities.
  Is = _indep_comp_cartesian_indices(SymTensorValue{D})
  ws = _stored_comp_weights(SymTensorValue{D})
  comps = [ _weighted_sum_expr(K -> :( a[$i,$j,$(K[1]),$(K[2])]*b[$(K[1]),$(K[2])] ), Is, ws)
            for (i,j) in Is ]
  :( SymTensorValue{Da}($(comps...)) )
end

# c_ij = a_kl*b_klij, the two free indices being the two last of `b`
@generated function double_contraction(a::AbstractSymTensorValue{Da,Ta}, b::SymFourthOrderTensorValue{Db,Tb}) where {Da,Db,Ta,Tb}
  can_contract(Da,Db) || return :( contracted_product(Val(2), a, b) ) # for error msg
  D = Int(Db)
  iszero(D) && return :(zero(SymTensorValue{Db,$(Base.promote_op(*,Ta,Tb))}))
  # `a` is contracted over the whole tensor, read off the stored comp indices
  # of the symmetric layout weighted by their multiplicities.
  Is = _indep_comp_cartesian_indices(SymTensorValue{D})
  ws = _stored_comp_weights(SymTensorValue{D})
  comps = [ _weighted_sum_expr(K -> :( a[$(K[1]),$(K[2])]*b[$(K[1]),$(K[2]),$i,$j] ), Is, ws)
            for (i,j) in Is ]
  :( SymTensorValue{Db}($(comps...)) )
end

# c_ij = a_ijkl*b_kl.
function double_contraction(a::SymFourthOrderTensorValue, b::MultiValue{Tuple{Db,Dc}}) where {Db,Dc}
  double_contraction(a, symmetric_part(b))
end

# c_ij = a_kl*b_klij
function double_contraction(a::MultiValue{Tuple{Da,Db}}, b::SymFourthOrderTensorValue) where {Da,Db}
  double_contraction(symmetric_part(a), b)
end

const ⋅² = double_contraction

###############################################################
# Reductions
###############################################################

for op in (:sum, :maximum, :minimum)
  @eval begin
    $op(a::MultiValue) = $op(a.data)
  end
end

###############################################################
# Outer product (aka tensor product, also dyadic product)
###############################################################

outer(a::Number, b::Number) = a * b

outer(a::MultiValue, b::Number) = a * b
outer(a::Number, b::MultiValue) = a * b

"""
    outer(a,b)
    a ⊗ b

Outer product (or tensor-product) of two `Number`s and/or `MultiValue`s, that is
`(a⊗b)[i₁,...,iₙ,j₁,...,jₙ] = a[i₁,...,iₙ]*b[j₁,...,jₙ]`. This falls back to standard
multiplication if `a` or `b` is a scalar.
"""
outer(a::MultiValue, b::MultiValue) = contracted_product(Val(0), a, b)

# c_ijkl = a_ij*b_kl
@generated function outer(a::AbstractSymTensorValue{D,Ta}, b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :(zero(SymFourthOrderTensorValue{D,$(Base.promote_op(*,Ta,Tb))}))
  comps = [ :( a[$i,$j]*b[$k,$l] )
            for (i,j,k,l) in _indep_comp_cartesian_indices(SymFourthOrderTensorValue{D}) ]
  :( SymFourthOrderTensorValue{D}($(comps...)) )
end

const ⊗ = outer

###############################################################
# Cross Product
###############################################################

function cross(a::MultiValue{Tuple{3}}, b::MultiValue{Tuple{3}})
  VectorValue{3}(a[2]b[3] - a[3]b[2], a[3]b[1] - a[1]b[3], a[1]b[2] - a[2]b[1])
end

function cross(a::MultiValue{Tuple{2}}, b::MultiValue{Tuple{2}})
  a[1]b[2] - a[2]b[1]
end

"""
    cross(a::VectorValue{3}, b::VectorValue{3}) -> VectorValue{3}
    cross(a::VectorValue{2}, b::VectorValue{2}) -> scalar
    a × b

Cross product of 2D and 3D vector.
"""
cross(a::MultiValue, b::MultiValue) = error("Cross product only defined for R2 and R3 vectors of same dimension")

###############################################################
# Linear Algebra
###############################################################

rank(a::MultiValue) = rank(SArray(a))

"""
    det(a::MultiValue{Tuple{D,D},T})

Determinent of square second order tensors.
"""
function det(a::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  D1 == D2 || @unreachable "det only defined for square tensors, got size $(size(a))"
  det(get_array(a))
end
det(a::MultiValue) = @unreachable "det only defined for square tensors, got size $(size(a))"

det(a::MultiValue{Tuple{1,1}}) = a[1]

function det(a::MultiValue{Tuple{2,2}})
  a_11 = a[1, 1]; a_12 = a[1, 2]
  a_21 = a[2, 1]; a_22 = a[2, 2]
  a_11*a_22 - a_12*a_21
end

function det(a::MultiValue{Tuple{3,3}})
  a_11 = a[1,1]; a_12 = a[1,2]; a_13 = a[1,3]
  a_21 = a[2,1]; a_22 = a[2,2]; a_23 = a[2,3]
  a_31 = a[3,1]; a_32 = a[3,2]; a_33 = a[3,3]
  a_11*a_22*a_33 + a_12*a_23*a_31 + a_13*a_21*a_32 -
    (a_11*a_23*a_32 + a_12*a_21*a_33 + a_13*a_22*a_31)
end

det(::SkewSymTensorValue{3,T}) where T = zero(T)

"""
    inv(a::MultiValue{Tuple{D1,D2}})

Inverse of a second order tensor.

`inv` inverse the variance of both indices.
"""
function inv(a::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  Int(D1) == Int(D2) || @unreachable "inv undefined for this tensor shape: $(size(a))"
  TensorValue{inverse_variance(D2),inverse_variance(D1)}(inv(get_array(a)))
end

"""
    pinvJt(Jt::MultiValue{Tuple{D1,D2}})

(Pseudo-)inverse of a transposed Jacobian `Jt`. This is the gradient
`iJt=∇φ⁻¹` of the inverse of a diffeomorphism `φ`, computed from its gradient
`Jt=∇φ`.

If `Jt` is full rank, that is `D = Int(min(D1,D2))` the manifold dimension,
`pinvJt(Jt)` is
- `inv(Jt)`     if `D1` = `D2`,
- `J⋅inv(Jt⋅J)` if `D1` < `D2`,
- `inv(J⋅Jt)⋅J` if `D1` > `D2`,
with `J = transpose(Jt)`. These closed forms are the Moore-Penrose
pseudo-inverse of `Jt`, to which `pinvJt` falls back when `Jt` is rank
deficient. Like the inverse, the pseudo-inverse swaps the position and
variance of both indices of the argument.

# Extended help

Let φ :  ℝᴰᵐ ⊃ M → N ⊂ ℝᴰⁿ be a ``C^1``-diffeomorphism between manifolds M
and N of dimension Dm = min(Dm,Dn), and ψ = φ⁻¹ its inverse. Write Jφ for the
Jacobian (transposed gradient) of φ, and Jψ for that of ψ. Both have a
contravariant first index, in the image of their map, and a covariant
second index, in its domain. Then

    Jφᵀ = pinvJt(Jψᵀ) = (Jψ⋅Jψᵀ)⁻¹⋅Jψ
    Jψᵀ = pinvJt(Jφᵀ) = Jφ⋅(Jφᵀ⋅Jφ)⁻¹

where these relations collapse to Jφ = Jψ⁻¹ if Dm == Dn. Remark that these
relations are not correct with regard to variance (two indices raising and
lowering are needed). This is because Jφ and Jψ are the matrices of the tangent
maps by φ and ψ and are thus coordinate dependent, they are not proper tensors.
Conversion so `SArray` are used for the computations.

Consider p̂ ∈ M, p = φ(p̂) ∈ N, and v̂ ∈ Tp̂M, ω̂ ∈ Tp̂*M, v ∈ TpN, ω ∈ Tp*M. The
four transports between the (co)tangent spaces of M and N are

    v = φ_* v̂ = ψ^* v̂,   vⁱ = v̂ʲ Jφ[j,i]
    ω̂ = ψ_* ω = φ^* ω,   ω̂ⱼ = Jφ[j,i] ωᵢ
    v̂ = ψ_* v = φ^* v,   v̂ʲ = vⁱ Jψ[i,j]
    ω = φ_* ω̂ = ψ^* ω̂,   ωᵢ = Jψ[i,j] ω̂ⱼ

The pushforward of a covariant tensor and the pullback of a contravariant one
are transports by the inverse map, hence the need for its inverse Jacobian and
`pinvJt`.
"""
function pinvJt(J_or_Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  @assert is_opposit_variance(D1,D2) "the argument is not the gradient or the
  Jacobian of a cell map, it needs be of type `(1,1)`, got
  $(tensor_type(J_or_Jt)[1:2])"

  V = TensorValue{inverse_variance(D2),inverse_variance(D1)}
  m = get_array(J_or_Jt)

  # pinv is not systematically used because it is very slow (500ns minimum, vs
  # 10-20 for the full rank formula) TODO: implement custom lu for small
  # MMatrix in StaticArrays to avoid calling blas, to optimize pinv.
  m_pinv = if D1 == D2
    iszero(abs(det(m))) ? pinv(m) : inv(m)
  elseif D1 < D2
    iszero(abs(det(m*m'))) ? pinv(m) : m'*inv(m*m')
  else # D1 > D2
    iszero(abs(det(m'*m))) ? pinv(m) : inv(m'*m)*m'
  end
  V(m_pinv)
end

const InverseStableTensorTypes{D} = Union{SymTensorValue{D},SkewSymTensorValue{D}}

for V in (:SymTensorValue, :SkewSymTensorValue)
  @eval function inv(a::$V{D}) where D
    ai = inv(get_array(a))
    $V{inverse_variance(D),eltype(ai)}(ai)
  end
end

function inv(a::MultiValue{Tuple{1,1}})
  r = 1 / a[1]
  T = change_eltype(a, typeof(r))
  T(r)
end

function inv(a::MultiValue{Tuple{2,2}})
  c = 1 / det(a)
  data = (a[2, 2] * c, -a[2, 1] * c, -a[1, 2] * c, a[1, 1] * c)
  TensorValue{2}(data)
end

function inv(a::MultiValue{Tuple{3,3}})
  a_11 = a[1,1]; a_12 = a[1,2]; a_13 = a[1,3]
  a_21 = a[2,1]; a_22 = a[2,2]; a_23 = a[2,3]
  a_31 = a[3,1]; a_32 = a[3,2]; a_33 = a[3,3]
  c = 1/det(a)
  data = (
     ( a_22*a_33 - a_23*a_32 )*c,
    -( a_21*a_33 - a_23*a_31 )*c,
     ( a_21*a_32 - a_22*a_31 )*c,
    -( a_12*a_33 - a_13*a_32 )*c,
     ( a_11*a_33 - a_13*a_31 )*c,
    -( a_11*a_32 - a_12*a_31 )*c,
     ( a_12*a_23 - a_13*a_22 )*c,
    -( a_11*a_23 - a_13*a_21 )*c,
     ( a_11*a_22 - a_12*a_21 )*c)
  TensorValue{3}(data)
end

function inv(a::SymTensorValue{2})
  c = 1/det(a)
  T = change_eltype(a,typeof(c))
  T(a[2,2]*c, -a[2,1]*c, a[1,1]*c)
end

inv(::SymTracelessTensorValue{1,T}) where T = TensorValue{1,1}(inv(zero(T)))
function inv(a::SymTracelessTensorValue{2})
  c = -1/det(a)
  T = change_eltype(a,typeof(c))
  T(a[1,1]*c, a[2,1]*c)
end

inv(::SkewSymTensorValue{1,T}) where T = TensorValue{1,1}(inv(zero(T)))
inv(a::SkewSymTensorValue{2}) = (typeof(a))(-inv(a.data[1]))
inv(a::SkewSymTensorValue{3,T,L}) where {T,L} = SkewSymTensorValue{3,T}(tfill(inv(zero(T)), Val(L)))


# TODO: write tests for the adjugate, based on adjugate(a)⋅a = det(a)*one(a), except for SymTracelessTensorValue that has no "one", but promotes to SymTensorValue.
"""
    adjugate(a::MultiValue{Tuple{D1,D2}})

The transpose of the co-factor tensor of `a`, such that

    adjugate(a)⋅a = det(a)*one(a).

It is well defined (no NaN) even if `det(a)` is zero.
`adjugate` inverses the variance of both indices.
"""
function adjugate(a::MultiValue{Tuple{D1,D2},T}) where {D1,D2,T}
  D1==D2 || throw(ArgumentError("adjugate is only defined on square tensors"))

  D1_inv, D2_inv = inverse_variance.((D1,D2))
  V = TensorValue{D1_inv,D2_inv}
  if     D1==0
    V{T}()
  elseif D1==1
    one(V{T})
  elseif D1==2
    V( a[2,2], -a[2,1],
      -a[1,2],  a[1,1])
  elseif D1==3
    c11 =  (a[2,2]*a[3,3] - a[2,3]*a[3,2])
    c12 = -(a[2,1]*a[3,3] - a[2,3]*a[3,1])
    c13 =  (a[2,1]*a[3,2] - a[2,2]*a[3,1])
    c21 = -(a[1,2]*a[3,3] - a[1,3]*a[3,2])
    c22 =  (a[1,1]*a[3,3] - a[1,3]*a[3,1])
    c23 = -(a[1,1]*a[3,2] - a[1,2]*a[3,1])
    c31 =  (a[1,2]*a[2,3] - a[1,3]*a[2,2])
    c32 = -(a[1,1]*a[2,3] - a[1,3]*a[2,1])
    c33 =  (a[1,1]*a[2,2] - a[1,2]*a[2,1])
    V( # this is transposed since Julia is column major
        c11, c21, c31,
        c12, c22, c32,
        c13, c23, c33
    )
  else
    _numerical_adjugate(V,a)
  end
end

function adjugate(a::SymTensorValue{D,T}) where {D,T}
  D_inv = inverse_variance(D)
  V = SymTensorValue{D_inv}
  if D==0
    V{T}()
  elseif D==1
    one(V{T})
  elseif D==2
    V(a[2,2], -a[2,1],
               a[1,1])
  elseif D==3
    c11 =  (a[2,2]*a[3,3] - a[2,3]^2)
    c12 = -(a[2,1]*a[3,3] - a[2,3]*a[3,1])
    c13 =  (a[2,1]*a[3,2] - a[2,2]*a[3,1])
    c22 =  (a[1,1]*a[3,3] - a[1,3]^2)
    c23 = -(a[1,1]*a[3,2] - a[1,2]*a[3,1])
    c33 =  (a[1,1]*a[2,2] - a[1,2]^2)
    return V(
        c11, c12, c13,
             c22, c23,
                  c33
    )
  else
    _numerical_adjugate(V,a)
  end
end

function adjugate(a::SymTracelessTensorValue{D,T}) where {D,T}
  adjugate(convert(SymTensorValue{D,T}, get_array(a)))
end

function _numerical_adjugate(::Type{V},a) where V
    adj_a = try inv(a)*det(a) # faster than the svd version, and non-allocating
    catch
      _svd_adjugate(a)
    end
    convert(V, adj_a)
end

function _svd_adjugate(a)
  # Credit: Steven G. Johnson (https://discourse.julialang.org/t/computing-the-adjugate-cofactors-of-a-matrix/125184)
  # This allocates a little bit
  function adjugate(D::Diagonal{<:Number})
    d = D.diag
    Base.require_one_based_indexing(d)
    length(d) < 2 && return Diagonal(one.(d))
    # compute diagonal of adj(D): dadj[i] = prod(d) / d[i], but avoid dividing by zero
    dadj = similar(d)
    prod = one(eltype(d))
    for i = 1:length(d)
      dadj[i] = prod
      prod *= d[i]
    end
    prod = one(eltype(d))
    for i = length(d):-1:1
      dadj[i] *= prod
      prod *= d[i]
    end
    return Diagonal(dadj)
  end

  F = svd(SMatrix(a))
  det_VtU = det(F.Vt * F.U)
  iszero(det_VtU) && return zero(a)
  F.V * (adjugate(Diagonal(F.S)) * det_VtU) * F.U'
end


"""
    eigen(a::MultiValue{Tuple{D1,D2}})

Eigenvalue decomposition of a square second order tensor.

Eigenvalues are coordinate-invariant only for a (1,1) tensor, so `D1` and `D2`
must be contractible, see [`can_contract`](@ref). The decomposition is returned
in terms of `SArray`s, which carry no variance.
"""
function eigen(a::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  can_contract(D1,D2) || @unreachable "eigen undefined for this tensor shape: $(size(a))"
  eigen(get_array(a))
end
eigen(a::MultiValue) = @unreachable "eigen undefined for this tensor shape: $(size(a))"

"""
    normalize(a::MultiValue)

Normalization of a tensor value.
"""
normalize(a::MultiValue) = a / norm(a)

"""
    sqrt(a::MultiValue{Tuple{D1,D2}})

Square root of a square second order tensor, carrying the indices of `a`.

A square root is coordinate-invariant only for a (1,1) tensor, so `D1` and `D2`
must be contractible, see [`can_contract`](@ref).
"""
function sqrt(t::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  can_contract(D1,D2) || @unreachable "sqrt undefined for this tensor shape: $(size(t))"
  TensorValue{D1,D2}(sqrt(SArray(t)))
end

###############################################################
# Measure
###############################################################

"""
    meas(a::MultiValue{Tuple{D}}, metric_args...)
    meas(a::MultiValue{Tuple{1,D2}})

Norm of a vector.
"""
meas(a::MultiValue{Tuple{D}}, metric_args::Vararg{_MetricArg}) where D =
  sqrt(inner(a, a, metric_args...))

function meas(v::MultiValue{Tuple{1,D}}) where D
  # The row is read as a plain vector, which is only a length-preserving
  # reading when `D` carries no variance.
  D isa Int || @unreachable "meas undefined for this tensor shape: $(size(v))"
  t = VectorValue(v.data)
  meas(t)
end

function meas(v::MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1, n2, n3)
  meas(n)
end

"""
    meas(J::MultiValue{Tuple{D1,D2}})

Returns the absolute `D1`-dimensional volume of the parallelepiped
formed by the rows of `J`, that is `sqrt(det(J⋅Jᵀ))`, or `abs(det(J))` if `D1`=`D2`.
This is used to compute the contribution of the metric `G = J⋅Jᵀ` of a changes of variables in integrals.

If it has variance, `J` may only be a cell map gradient of size `(Cov(D1), Con(D2))`.
"""
@generated function meas(Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  if D1 isa Cov && D2 isa Con
    V = TensorValue{Int(D1),Int(D2)}
    return :( meas($V(get_array(Jt))) )
  end

  if !(D1 isa Int && D2 isa Int)
    msg = "meas undefined for this tensor shape: $(size(Jt)), it is not the
      transposed Jacobian of a chart"
    return :( @unreachable $msg )
  end

  D1 == D2 && return :( abs(det(Jt)) )
  :( sqrt(det(Jt ⋅ transpose(Jt))) )
end

"""
    norm(a::MultiValue, metric_args...)
    norm(a::MultiValue, p::Number, metric_args...)

2-norm of `a` (`p` must be 2), defined as `sqrt(inner(u,u,metric_args...))`.

See also [`parse_metric_data`](@ref) and [`inner`](@ref).
"""
@inline norm(u::MultiValue{Tuple{D},<:Real}, metric_args::Vararg{_MetricArgNoDetG}) where D =
  sqrt(inner(u, u, metric_args...))
@inline norm(u::MultiValue{Tuple{D}}, metric_args::Vararg{_MetricArgNoDetG}) where D =
  sqrt(real(inner(u, conj(u), metric_args...)))
@inline norm(u::MultiValue{Tuple{D1,D2},<:Real}, metric_args::Vararg{_MetricArgNoDetG}) where {D1,D2} =
  sqrt(inner(u, u, metric_args...))
@inline norm(u::MultiValue{Tuple{D1,D2}}, metric_args::Vararg{_MetricArgNoDetG}) where {D1,D2} =
  sqrt(real(inner(u, conj(u), metric_args...)))
@inline norm(u::MultiValue{Tuple{0},T}) where T<:Real= sqrt(zero(T))
@inline norm(u::MultiValue{Tuple{0},T}) where T = sqrt(real(zero(T)))

@inline function norm(u::MultiValue, p::Real, metric_args::Vararg{_MetricArgNoDetG})
  p == 2 || throw(ArgumentError("only the 2-norm of a MultiValue is defined, got p = $p"))
  norm(u, metric_args...)
end

###############################################################
# conj, real, imag
###############################################################

for op in (:conj,:real,:imag)
  @eval begin
    function ($op)(a::T) where T<:MultiValue
      Li = num_indep_components(a)
      r = map($op, Tuple(a)[1:Li])
      T2 = _eltype($op,r,a)
      M = change_eltype(a,T2)
      M(r)
    end
  end
end

###############################################################
# Trace
###############################################################

# Rejection of a pair of indices that cannot be traced over.
_not_traceable_msg(i, j) = "Second order tensor of indices $i and $j is not square, see `can_contract`."
_not_traceable_expr(i, j) = :( throw(ArgumentError($(_not_traceable_msg(i,j)))) )

"""
    tr(v::MultiValue{Tuple{D1,D2}})

Return the trace of a second order square tensor, defined by `Σᵢ vᵢᵢ` or 0 if `D1`=0.

The traced indices must be contractible, see [`can_contract`](@ref).
"""
@generated function tr(v::MultiValue{Tuple{D1,D2},T}) where {D1,D2,T}
  can_contract(D1,D2) || return _not_traceable_expr(D1,D2)
  D = Int(D1)
  iszero(D) && return :(zero(T))
  summands = [:(v[$i,$i]) for i in 1:D]
  Expr(:call, :+, summands...)
end

for V in (:SymTracelessTensorValue, :SkewSymTensorValue)
  @eval @generated function tr(v::$V{D,T}) where {D,T}
    can_contract(D,D) || return _not_traceable_expr(D,D)
    :(zero(T))
  end
end

"""
    tr(v::MultiValue{Tuple{D,D,D2}}) -> ::VectorValue{D2}

Return a vector of length `D2` of traces computed on the first two indices: `resⱼ = Σᵢ vᵢᵢⱼ`.
"""
@generated function tr(v::NTensorValue{S,T,N}) where {S,T,N}
  if N > 3
    # This may be generalizable if wanted
    msg = "Trace is not implemented for tensors of order ≥ 4, got tensors ::$(v)"
    return :( throw(ArgumentError($msg)) )
  end

  A12 = first(size(v), 2)
  can_contract(A12...) || return _not_traceable_expr(A12...)

  B = last(size(v))
  V = _tensor_base_type_expr(1, (B,))
  iszero(length(v)) && return :(zero($V{T}))

  comps = Expr[]
  for k in 1:Int(B)
    summands = [:( v[$i,$i,$k] ) for i in 1:Int(A12[1])]
    push!(comps, Expr(:call, :+, summands...))
  end
  :( $V($(comps...)) )
end

###############################################################
# Adjoint and transpose
###############################################################

adjoint(a::MultiValue{Tuple{D,D}}) where D = @notimplemented
transpose(a::MultiValue{Tuple{D,D}}) where D = @notimplemented

@generated function adjoint(a::TensorValue{D1,D2,T}) where {D1,D2,T}
  comps = Expr[]
  for i in 1:Int(D1)
    for j in 1:Int(D2)
      push!(comps, :(conj(a[$i,$j])) )
    end
  end
  :( TensorValue{D2,D1,T}($(comps...)) )
end

@generated function transpose(a::TensorValue{D1,D2,T}) where {D1,D2,T}
  comps = Expr[]
  for i in 1:Int(D1)
    for j in 1:Int(D2)
      push!(comps, :(a[$i,$j]) )
    end
  end
  :( TensorValue{D2,D1,T}($(comps...)) )
end

@inline function adjoint(a::TensorValue{D1,D2,T}) where {D1,D2,T<:Real}
  transpose(a)
end

adjoint(a::AbstractSymTensorValue) = conj(a)
adjoint(a::SkewSymTensorValue) = -conj(a)

transpose(a::AbstractSymTensorValue) = a
transpose(a::SkewSymTensorValue) = -a

###############################################################
# permutedims
###############################################################

"""
    permutedims(a::MultiValue{Sa,Ta,Na}, perm::NTuple{Na,Int})
    permutedims(a::MultiValue{Sa,Ta,Na},     ::Val{perm})

Return a tensor whose indices are those of `a` reordered by `perm`:
`result[i₁,…,iₙ] = a[i_{σ⁻¹(1)},…,i_{σ⁻¹(n)}]`, where `σ = perm`.

!!! warning
    The method accepting a `Tuple` is not type inferable and it allocates, pass
    the `perm` by `Val` if it is known at compile time.

The output shape is `(size(a, perm[1]), …, size(a, perm[N]))`.
For second-order tensors this is equivalent to `transpose` (with `perm = (2,1)`).
Symmetry of the input tensor is not preserved in the output type.

```@example
t = TensorValue{2,3}(1:6...)
permutedims(t, Val((2,1)) ) # 2.123 ns  (0 allocations: 0 bytes)
permutedims(t, (2,1))       # 82.367 ns (3 allocations: 160 bytes)
# -> TensorValue{3, 2, Int64, 6}(1, 3, 5, 2, 4, 6)
```
"""
function Base.permutedims(
  a::MultiValue{Sa,Ta,Na},
  perm::NTuple{Na,Int}) where {Sa,Ta,Na}
  permutedims(a, Val(perm))
end

@generated function Base.permutedims(
  a::MultiValue{Sa,Ta,Na}, ::Val{P}) where {Sa,Ta,Na,P}

  @assert P isa Tuple || P isa Int

  sort([P...]) == collect(1:Na) || begin
    msg = "$P is not a valid permutation of 1:$Na"
    return :(throw(ArgumentError($msg)))
  end

  sz_a = size(a)
  Sr = ntuple(k -> sz_a[P[k]], Na)
  V = _tensor_base_type_expr(Na,Sr)

  iszero(length(a)) && return :( zero($V) )

  inv_perm = zeros(Int, Na)
  for k in 1:Na
    inv_perm[P[k]] = k
  end

  comps = Expr[]
  for ci in CartesianIndices(Int.(Sr))
    a_idx = ntuple(k -> ci[inv_perm[k]], Na)
    push!(comps, :( a[$a_idx...] ) )
  end
  :( $V($(comps...)) )
end

###############################################################
# Symmetric and Skew symmetric parts
###############################################################

"""
    symmetric_part(v::MultiValue{Tuple{D,D}})::AbstractSymTensorValue

Return the symmetric part of second order tensor, that is `½(v + vᵀ)`.
Return `v` if  `v isa AbstractSymTensorValue`, and the zero symmetric tensor if
`v  isa SkewSymTensorValue`.
"""
@generated function symmetric_part(v::MultiValue{Tuple{D1,D2},T}) where {D1,D2,T}
  D1 === D2 || return :( @notimplemented $(
    "symmetric_part requires two identical indices, got $D1 and $D2") )
  D = Int(D1)
  iszero(D) && return :(zero(SymTensorValue{D1,T}))
  comps = Expr[]
  for j in 1:D
    for i in j:D
      push!(comps, :((v[$i,$j] + v[$j,$i])/2) )
    end
  end
  :( SymTensorValue{D1}($(comps...)) )
end

symmetric_part(v::AbstractSymTensorValue) = v
symmetric_part(::SkewSymTensorValue{D,T}) where {D,T} = zero(SymTensorValue{D,T})

"""
    skew_symmetric_part(v::MultiValue{Tuple{D,D}})::SkewSymTensorValue{D}

Return the asymmetric part of `v`, that is `½(v - vᵀ)`.
Return the zero skew symmetric tensor if `v  isa AbstractSymTensorValue`, and
`v` itself if  `v  isa SkewSymTensorValue`.
"""
@generated function skew_symmetric_part(v::MultiValue{Tuple{D1,D2},T}) where {D1,D2,T}
  D1 === D2 || return :( @notimplemented $(
    "skew_symmetric_part requires two identical indices, got $D1 and $D2") )
  D = Int(D1)
  iszero(D) && return :(zero(SkewSymTensorValue{D1,T}))
  comps = Expr[]
  for i in 1:D
    for j in i+1:D
      push!(comps, :((v[$i,$j] - v[$j,$i])/2) )
    end
  end
  :( SkewSymTensorValue{D1}($(comps...)) )
end

skew_symmetric_part(::AbstractSymTensorValue{D,T}) where {D,T} = zero(SkewSymTensorValue{D,T})
skew_symmetric_part(v::SkewSymTensorValue) = v

###############################################################
# diag
###############################################################

function LinearAlgebra.diag(a::MultiValue{Tuple{D,D},T}) where {D,T}
  D isa Int || @unreachable "diag requires un-varianced indices, got $D"
  VectorValue{D,T}((a[i,i] for i in 1:D)...)
end

function LinearAlgebra.diag(a::SkewSymTensorValue{D,T}) where {D,T}
  D isa Int || @unreachable "diag requires un-varianced indices, got $D"
  zero(VectorValue{D,T})
end

###############################################################
# Broadcast
###############################################################

# Tensor components need to be relative to comparable tensor bases for
# operations to be well defines.
function _broadcast_check(a::MultiValue, b::MultiValue)
  size(a) === size(b) || throw(DimensionMismatch(
    "Cannot broadcast over tensors of sizes $(size(a)) and $(size(b))"))
end

function Base.broadcasted(f, a::VectorValue, b::VectorValue)
  _broadcast_check(a, b)
  VectorValue(map(f, a.data, b.data))
end

function Base.broadcasted(f, a::ContraVectorValue{D}, b::ContraVectorValue) where D
  _broadcast_check(a, b)
  ContraVectorValue{D}(map(f, a.data, b.data))
end

function Base.broadcasted(f, a::TensorValue{D1,D2}, b::TensorValue)  where {D1,D2}
  _broadcast_check(a, b)
  TensorValue{D1,D2}(map(f, a.data, b.data))
end

function Base.broadcasted(f, a::AbstractSymTensorValue{D}, b::AbstractSymTensorValue) where D
  _broadcast_check(a, b)
  SymTensorValue{D}(map(f, a.data, b.data))
end
