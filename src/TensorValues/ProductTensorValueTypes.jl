###############################################################
# ProductTensorValue helpers
###############################################################

"""
    _minimal_tensortype_with_shape(V::Type{<:MultiValue{S}})

Return `V` with all type parameter removed, except those necessary for the S
parameter of the MultiValue supertype to still be defined.
"""
_minimal_tensortype_with_shape(::Type{<:ScalarValue})                                = ScalarValue
_minimal_tensortype_with_shape(::Type{<:VectorValue{D}})              where D        = VectorValue{D}
_minimal_tensortype_with_shape(::Type{<:TensorValue{D}})              where D        = TensorValue{D,D}
_minimal_tensortype_with_shape(::Type{<:TensorValue{D1,D2}})          where {D1,D2}  = TensorValue{D1,D2}
_minimal_tensortype_with_shape(::Type{<:SymTensorValue{D}})           where D        = SymTensorValue{D}
_minimal_tensortype_with_shape(::Type{<:SymTracelessTensorValue{D}})  where D        = SymTracelessTensorValue{D}
_minimal_tensortype_with_shape(::Type{<:SkewSymTensorValue{D}})       where D        = SkewSymTensorValue{D}
_minimal_tensortype_with_shape(::Type{<:SymFourthOrderTensorValue{D}})where D        = SymFourthOrderTensorValue{D}
_minimal_tensortype_with_shape(::Type{<:ThirdOrderTensorValue{D}})    where D        = ThirdOrderTensorValue{D,D,D}
_minimal_tensortype_with_shape(::Type{<:ThirdOrderTensorValue{D1,D2,D3}}) where {D1,D2,D3} = ThirdOrderTensorValue{D1,D2,D3}
_minimal_tensortype_with_shape(::Type{<:NTensorValue{S}})     where S        = NTensorValue{S}
_minimal_tensortype_with_shape(::Type{<:SymNTensorValue{S}})          where S        = SymNTensorValue{S}
_minimal_tensortype_with_shape(::Type{<:SkewSymKTensorValue{V,K,D}})  where {V,K,D}  = SkewSymKTensorValue{V,K,D}
_minimal_tensortype_with_shape(::Type{V}) where V<:MultiValue = @unreachable "$V cannot be a
  factor of a ProductTensorValue, `_minimal_tensortype_with_shape` is not implemented for it."

# Stored and independent component counts of both factors of the product of `A`
# and `B`: `(La, Lb, iLa, iLb)`. `A` and `B` must be canonical factor argument
# types, so that these counts describe the same tensors as the product size.
function _prod_counts(::Type{A}, ::Type{B}) where {A,B}
  La  = num_stored_components(A)
  Lb  = num_stored_components(B)
  iLa = num_indep_components(A)
  iLb = num_indep_components(B)
  (La, Lb, iLa, iLb)
end

# Canonical factor argument types, size, order, and number of stored and of independent
# components of the product of `A` and `B`.
function _prod_shape_impl(::Type{A}, ::Type{B}) where {A,B}
  Ac, Bc = _minimal_tensortype_with_shape(A), _minimal_tensortype_with_shape(B)
  S = Tuple{size(Ac)..., size(Bc)...}
  La, Lb, iLa, iLb = _prod_counts(Ac,Bc)
  N = length(_Size_to_size(S))
  (Ac, Bc, S, N, La*Lb, iLa*iLb)
end

# @generated ensures compile time constants and inference of the parameters
@generated _prod_shape(::Type{A}, ::Type{B}) where {A,B} = :( $(_prod_shape_impl(A,B)) )

# The stored components of a product, from its independent ones.
@generated function _prod_stored_components(
  ::Type{A}, ::Type{B}, ::Type{T}, data::Tuple) where {A,B,T}

  La, Lb, iLa, iLb = _prod_counts(A,B)
  La*Lb == iLa*iLb && return :( data )

  ca = Tuple.(component_basis(make_concretetype(A{Int})))
  cb = Tuple.(component_basis(make_concretetype(B{Int})))

  comps = Expr[]
  for kb in 1:Lb, ka in 1:La
    terms = Expr[]
    for lb in 1:iLb, la in 1:iLa
      c = ca[la][ka]*cb[lb][kb]
      iszero(c) && continue
      l = (lb-1)*iLa + la
      push!(terms, isone(c)  ? :(  data[$l] ) :
                   isone(-c) ? :( -data[$l] ) : :( $c*data[$l] ))
    end
    push!(comps, isempty(terms) ? :( zero(T) ) : Expr(:call, :+, terms...))
  end
  Expr(:tuple, comps...)
end

"""
    _indep_comp_cartesian_indices(::Type{V<:MultiValue})

A full tensor index at which each independent component of `V` can be read, that
is the position where `V`'s component basis reduces to the `l`th Cartesian basis
vector, so that for any `u::V`

    get_array(u)[_indep_comp_cartesian_indices(V)[l]...] == indep_comp_getindex(u,l)
"""
function _indep_comp_cartesian_indices(::Type{V}) where V<:MultiValue
  VV = make_concretetype(_minimal_tensortype_with_shape(V){Int})
  basis = map(get_array, component_basis(VV))
  map(eachindex(basis)) do l
    for I in CartesianIndices(first(basis))
      all(basis[i][I] == (i == l ? 1 : 0) for i in eachindex(basis)) && return Tuple(I)
    end
    @unreachable "No component of $V holds its independent component $l alone, it
      cannot be read from a full tensor."
  end
end

"""
    _stored_comp_cartesian_indices(::Type{V<:MultiValue})

A full tensor index at which each *stored* component of `V` can be read, in the
order they are stored, so that for any `u::V`

    get_array(u)[_stored_comp_cartesian_indices(V)[k]...] == Tuple(u)[k]

A type that stores only independent components reads them at
[`_indep_comp_cartesian_indices`](@ref). The general case inverts
[`stored_comp_indices`](@ref), which holds `+k` exactly at the indices reading
the stored component `k`.
"""
function _stored_comp_cartesian_indices(::Type{V}) where V<:MultiValue
  V_sci = get_array(stored_comp_indices(V))
  map(1:num_stored_components(V)) do k
    for I in CartesianIndices(V_sci)
      V_sci[I] == k && return Tuple(I)
    end
    @unreachable "The $k th stored component of $V is not a component of the
    full tensor."
  end
end

# The component type given by a factor type, `nothing` when it leaves it free.
_factor_eltype(::Type{<:MultiValue{S,T}}) where {S,T} = T
_factor_eltype(::Type{<:MultiValue}) = nothing
_factor_eltype(::Type) = nothing

# The component type that the factor types impose on their product, `nothing`
# when neither gives one.
function _prod_eltype(::Type{A}, ::Type{B}) where {A,B}
  Ta, Tb = _factor_eltype(A), _factor_eltype(B)
  isnothing(Ta) && return Tb
  isnothing(Tb) && return Ta
  Base.promote_op(*, Ta, Tb)
end
_prod_eltype(::Type{A}, ::Type{B}, ::Type{Td}) where {A,B,Td} = something(_prod_eltype(A,B), Td)


###############################################################
# ProductTensorValue type
###############################################################

"""
    ProductTensorValue{A,B,T,S,N,L} <: MultiValue{S,T,N,L}

Tensor type used to store the tensor (outer) product of `a::A` and `b::B`,
without more redundancy than `a` and `b`'s storage (when they have dependent
components).

- `S` concatenates `Sa` and `Sb`,
- `N` is the sum `Na + Nb`
- `L` is the product `La`×`Lb`, where does are the numer of stored components

The product tensor has `num_indep_components(b)`×`num_indep_components(b)`
independent components.

The independent components of `A` vary fastest in the internal storage,
consistently with the column major ordering of Julia:

    p = product_tensor(a,b)
    La = num_indep_components(a)
    indep_comp_getindex(p, (lb-1)*La + la) == indep_comp_getindex(a,la)*indep_comp_getindex(b,lb)

"""
struct ProductTensorValue{A,B,T,S,N,L} <: MultiValue{S,T,N,L}
  data::NTuple{L,T}

  function ProductTensorValue{A,B}(data::NTuple{Lid,Td}) where {A,B,Lid,Td}
    Ac, Bc, S, N, L, Li = _prod_shape(A,B)
    @check Lid == Li "Invalid number of components in ProductTensorValue{$A,$B} constructor, expecting $Li, got $Lid"
    T = _prod_eltype(A,B,Td)
    new{Ac,Bc,T,S,N,L}(_prod_stored_components(Ac,Bc,T,NTuple{Li,T}(data)))
  end

  function ProductTensorValue{A,B,T,S,N,L}(data::Tuple) where {A,B,T,S,N,L}
    Ac, Bc, Sc, Nc, Lc, Li = _prod_shape(A,B)
    @check (S,N,L) === (Sc,Nc,Lc) "The product of $A and $B has size $Sc, order
      $Nc and $Lc stored components, got $S, $N and $L"
    @check length(data) == Li "Invalid number of components in ProductTensorValue{$A,$B} constructor, expecting $Li, got $(length(data))"
    Tt = _prod_eltype(A,B)
    @check isnothing(Tt) || promote_type(Tt,T) === T "The components of the
      product of $A and $B have type $Tt, which does not promote to $T"
    new{Ac,Bc,T,Sc,Nc,Lc}(_prod_stored_components(Ac,Bc,T,NTuple{Li,T}(data)))
  end
end

# prevent nested ProductTensorValue
_minimal_tensortype_with_shape(::Type{V}) where V<:ProductTensorValue = @unreachable "$V cannot be
  a factor of a ProductTensorValue, products do not nest. Use `outer` to combine
  more than two tensors."

###############################################################
# Constructors (ProductTensorValue)
###############################################################

ProductTensorValue{A,B}(data::Tuple) where {A,B} = ProductTensorValue{A,B}(promote(data...))

# An empty tuple carries no component type, only the factors can give it. The
# shape is read first, so that an invalid factor is reported as such.
function ProductTensorValue{A,B}(data::Tuple{}) where {A,B}
  _, _, S, N, L = _prod_shape(A,B)
  T = _prod_eltype(A,B)
  isnothing(T) && @unreachable "The component type of an empty
    ProductTensorValue{$A,$B} is given neither by its (empty) data nor by its
    factors, use the ProductTensorValue{A,B,T} constructor."
  ProductTensorValue{A,B,T,S,N,L}(data)
end
function ProductTensorValue{A,B,T}(data::Tuple) where {A,B,T}
  _, _, S, N, L = _prod_shape(A,B)
  ProductTensorValue{A,B,T,S,N,L}(data)
end
function ProductTensorValue{A,B,T,S}(data::Tuple) where {A,B,T,S}
  _, _, _, N, L = _prod_shape(A,B)
  ProductTensorValue{A,B,T,S,N,L}(data)
end
function ProductTensorValue{A,B,T,S,N}(data::Tuple) where {A,B,T,S,N}
  _, _, _, _, L = _prod_shape(A,B)
  ProductTensorValue{A,B,T,S,N,L}(data)
end

# ProductTensorValue Vararg constructor

ProductTensorValue{A,B}(data::Number...) where {A,B} = ProductTensorValue{A,B}(data)
ProductTensorValue{A,B,T}(data::Number...) where {A,B,T} = ProductTensorValue{A,B,T}(data)
ProductTensorValue{A,B,T,S}(data::Number...) where {A,B,T,S} = ProductTensorValue{A,B,T,S}(data)
ProductTensorValue{A,B,T,S,N}(data::Number...) where {A,B,T,S,N} = ProductTensorValue{A,B,T,S,N}(data)
ProductTensorValue{A,B,T,S,N,L}(data::Number...) where {A,B,T,S,N,L} = ProductTensorValue{A,B,T,S,N,L}(data)

# single AbstractArray argument constructors

# Independent components of a product, read from its full tensor, without
# validating unicity of values of an assumed independent component.
@generated function _flatten_prod_components(
  data::AbstractArray, ::Type{A}, ::Type{B}) where {A,B}

  Ac, Bc, = _prod_shape(A,B)
  s = Int.((size(Ac)..., size(Bc)...))
  msg = "expected an array of size $s"
  Ia, Ib = _indep_comp_cartesian_indices(Ac), _indep_comp_cartesian_indices(Bc)
  comps = [ :( data[$(cia...), $(cib...)] ) for cib in Ib for cia in Ia ]
  quote
    Base.require_one_based_indexing(data)
    @check $s == size(data) $msg
    return $(Expr(:tuple, comps...))
  end
end

ProductTensorValue{A,B}(data::AbstractArray{T}) where {A,B,T} =
  ProductTensorValue{A,B,_prod_eltype(A,B,T)}(_flatten_prod_components(data,A,B))
ProductTensorValue{A,B,T}(data::AbstractArray) where {A,B,T} =
  ProductTensorValue{A,B,T}(_flatten_prod_components(data,A,B))
ProductTensorValue{A,B,T,S}(data::AbstractArray) where {A,B,T,S} =
  ProductTensorValue{A,B,T,S}(_flatten_prod_components(data,A,B))
ProductTensorValue{A,B,T,S,N}(data::AbstractArray) where {A,B,T,S,N} =
  ProductTensorValue{A,B,T,S,N}(_flatten_prod_components(data,A,B))
ProductTensorValue{A,B,T,S,N,L}(data::AbstractArray) where {A,B,T,S,N,L} =
  ProductTensorValue{A,B,T,S,N,L}(_flatten_prod_components(data,A,B))

"""
    product_tensor(a::Number, b::Number)

The tensor product of `a` and `b`, like `outer(a,b)`, but returns a
[`ProductTensorValue`](@ref) preserving the dependent components information
when both factors are tensors.
"""
@generated function product_tensor(a::MultiValue, b::MultiValue)
  A, B, S, N, L, Li = _prod_shape(a,b)
  if iszero(Li)
    T = _prod_eltype(a,b)
    return :( ProductTensorValue{$A,$B,$T,$S,$N,$L}(()) )
  end

  iLa = num_indep_components(A)
  comps = [ :( indep_comp_getindex(a,$la) * indep_comp_getindex(b,$lb) )
            for lb in 1:(Li÷iLa) for la in 1:iLa ]
  :( ProductTensorValue{$A,$B}($(Expr(:tuple, comps...))) )
end

product_tensor(a::Number, b::Number) = a*b

###############################################################
# Conversions (ProductTensorValue)
###############################################################

@generated function _ProductTensorValue_to_array(arg::ProductTensorValue{A,B,T}) where {A,B,T}
  s = Int.(size(arg))
  comps = [ :( arg[$(Tuple(ci)...)] ) for ci in CartesianIndices(s) ]
  :( return SArray{$(Tuple{s...}),T}($(Expr(:tuple, comps...))) )
end

# Inverse conversion
convert(::Type{<:SArray{S2,T2}}, arg::ProductTensorValue) where {S2,T2} = SArray{S2,T2}(_ProductTensorValue_to_array(arg))
convert(::Type{<:MArray{S2,T2}}, arg::ProductTensorValue) where {S2,T2} = MArray{S2,T2}(_ProductTensorValue_to_array(arg))

# Internal conversion
convert(::Type{<:ProductTensorValue{A,B,T,S,N,L}}, arg::ProductTensorValue{A,B,T2,S,N,L}) where {A,B,T,T2,S,N,L} =
  ProductTensorValue{A,B,T,S,N,L}(get_indep_components(arg))
convert(::Type{<:ProductTensorValue{A,B,T,S,N,L}}, arg::ProductTensorValue{A,B,T,S,N,L}) where {A,B,T,S,N,L} = arg

# Conversions from and to the tensors holding the same components
#
# It is assumed that all full rank tensors of the TensorValues module store
# their components "column" major, so that two full rank tensors of the same
# size hold the same components in the same order

"""
    is_full_rank(V::Type{<:Number})
    is_full_rank(a::Number)

Wether all components of `a::V` are independent.
"""
is_full_rank(::Type{V}) where V<:MultiValue = length(V) == num_indep_components(V)
is_full_rank(::Type{<:Number}) = true
is_full_rank(a::Number) = is_full_rank(typeof(a))

# From ProductTensorValue

function convert(
  ::Type{V}, arg::ProductTensorValue{A,B,T2,S}) where {A,B,T2,S,T,V<:MultiValue{S,T}}

  is_full_rank(V) && is_full_rank(arg) ||
    (A === ScalarValue && V<:B) || (B === ScalarValue && V<:A) ||
    throw(MethodError(convert, (V, arg)))
  V(get_indep_components(arg))
end

# To ProductTensorValue

function _convert_to_product(
  PT::Type{<:ProductTensorValue{A,B,T}}, arg::MultiValue{Sa}) where {A,B,T,Sa}

  Va = typeof(arg)
  _, _, S = _prod_shape(A,B)
  Sa === S && (is_full_rank(PT) && is_full_rank(arg) ||
    (A === ScalarValue && Va<:B) || (B === ScalarValue && Va<:A)) ||
    throw(MethodError(convert, (PT, arg)))
  PT(get_indep_components(arg))
end

convert(PT::Type{<:ProductTensorValue{A,B,T}}, arg::MultiValue) where {A,B,T} =
  _convert_to_product(PT,arg)

convert(
  PT::Type{<:ProductTensorValue{A1,B1,T1,S}},
  arg::ProductTensorValue{A2,B2,T2,S}) where {A1,B1,T1,A2,B2,T2,S} =
  _convert_to_product(PT,arg)

# Promotion with full rank tensors

function promote_rule(
  ::Type{<:ProductTensorValue{A,B,Ta,S,N,L}},
  ::Type{<:ProductTensorValue{A,B,Tb,S,N,L}}) where {A,B,Ta,Tb,S,N,L}

  ProductTensorValue{A,B,promote_type(Ta,Tb),S,N,L}
end

function promote_rule(
  ::Type{<:ProductTensorValue{A,B,Ta,S,N,L}},
  ::Type{V}) where {A,B,Ta,S,N,L,Tb,V<:MultiValue{S,Tb}}

  V<:ProductTensorValue && return Union{} # V cannot be PTV{A,B,Tb}, this dispatches above
  is_full_rank(A) && is_full_rank(B) && is_full_rank(V) || return Union{}
  ProductTensorValue{A,B,promote_type(Ta,Tb),S,N,L}
end

###############################################################
# Comparison (ProductTensorValue)
###############################################################

(==)(a::ProductTensorValue{A,B,Ta,S}, b::ProductTensorValue{A,B,Tb,S}) where {A,B,Ta,Tb,S} = a.data == b.data
(==)(a::ProductTensorValue{Aa,Ba,Ta,S}, b::ProductTensorValue{Ab,Bb,Tb,S}) where {Aa,Ba,Ta,Ab,Bb,Tb,S} = get_array(a) == get_array(b)
(==)(a::ProductTensorValue{A,B,T,S}, b::MultiValue{S}) where {A,B,T,S} = get_array(a) == get_array(b)
(==)(a::MultiValue{S}, b::ProductTensorValue{A,B,T,S}) where {A,B,T,S} = get_array(a) == get_array(b)

###############################################################
# Other constructors and conversions (ProductTensorValue)
###############################################################

change_eltype(::Type{<:ProductTensorValue{A,B,T1}},::Type{T2}) where {A,B,T1,T2} =
  ProductTensorValue{A,B,T2}
change_eltype(::Type{<:ProductTensorValue{A,B,T1,S,N,L}},::Type{T2}) where {A,B,T1,T2,S,N,L} =
  ProductTensorValue{A,B,T2,S,N,L}

function change_variance(::Type{<:ProductTensorValue{A,B}}, ::Val{Sv}) where {A,B,Sv}
  _check_change_variance(ProductTensorValue{A,B}, Sv)
  Na = length(size(A))
  ProductTensorValue{change_variance(A, Val(Sv[1:Na])),
                     change_variance(B, Val(Sv[Na+1:end]))}
end

change_variance(::Type{<:ProductTensorValue{A,B,T}}, ::Val{Sv}) where {A,B,T,Sv} =
  (change_variance(ProductTensorValue{A,B}, Val(Sv))){T}

size(::Type{<:ProductTensorValue{A,B}}) where {A,B} =
  tuple(size(A)..., size(B)...)

length(::Type{<:ProductTensorValue{A,B}}) where {A,B} = length(A)*length(B)

num_indep_components(::Type{<:ProductTensorValue{A,B}}) where {A,B} =
  num_indep_components(A)*num_indep_components(B)

num_stored_components(::Type{<:ProductTensorValue{A,B}}) where {A,B} =
  num_stored_components(A)*num_stored_components(B)

function zero(::Type{<:ProductTensorValue{A,B,T}}) where {A,B,T}
  Ac, Bc, S, N, L, Li = _prod_shape(A,B)
  ProductTensorValue{Ac,Bc,T,S,N,L}(ntuple(_ -> zero(T), Val(Li)))
end

###############################################################
# Indexing helpers (ProductTensorValue)
###############################################################

# The multiplicity of each stored component of `V` in its full tensor, that is
# the number of full tensor indices that read it.
function _stored_comp_weights(V::Type{<:MultiValue})
  V_sci = stored_comp_indices(V)
  AV_sci = get_array(V_sci)
  map(k -> count(I -> abs(AV_sci[I]) == k, CartesianIndices(AV_sci)),
      1:num_stored_components(V))
end

@propagate_inbounds function _get_data(p::ProductTensorValue{A,B}, i::Integer) where {A,B}
  La, _, iLa, _ = _prod_counts(A,B)
  La == iLa && return p.data[i]
  lb, la = divrem(i-1, iLa)
  p.data[lb*La + la + 1]
end

@generated function get_indep_components(p::ProductTensorValue{A,B}) where {A,B}
  La, _, iLa, iLb = _prod_counts(A,B)
  Expr(:tuple, ( :( p.data[$((lb-1)*La + la)] ) for lb in 1:iLb for la in 1:iLa )...)
end

###############################################################
# VTK export (ProductTensorValue)
###############################################################

function indep_components_names(::Type{<:ProductTensorValue{A,B}}) where {A,B}
  names_a, names_b = indep_components_names(A), indep_components_names(B)
  iLa, iLb = num_indep_components(A), num_indep_components(B)

  if length(names_a) != iLa || length(names_b) != iLb
    @warn "indep_components_names is inconsistent with num_indep_components for
    a factor of ProductTensorValue{$A,$B}, falling back to default names [string(i) for i in 1:L]"
    return [string(i) for i in 1:(iLa*iLb)]
  end

  [ na*nb for nb in names_b for na in names_a ]
end

###############################################################
# Factor slicing (ProductTensorValue)
###############################################################

_first_factor(::Type{<:ProductTensorValue{A,B}}) where {A,B} = A
_first_factor(a::ProductTensorValue) = _first_factor(typeof(a))
_second_factor(::Type{<:ProductTensorValue{A,B}}) where {A,B} = B
_second_factor(a::ProductTensorValue) = _second_factor(typeof(a))

# `p` is the sum of `aₗ ⊗ bₗ` over the independent components of one factor, so
# slicing it along either factor gives the other one's, and an operation acting
# on one index block alone acts on those slices.
#
# A factor with no component leaves the type of the slices undetermined, as an
# empty tuple of them has no element type.
_no_factor_slice_msg(A, B) = :( @notimplemented "The factors of a
  ProductTensorValue are only accessible when both of them have at least one
  component, got $($A) and $($B)" )

# The second factors of `p`, one per independent component of its first factor.
@generated function _second_factors(p::ProductTensorValue{A,B,T}) where {A,B,T}
  _, _, iLa, iLb = _prod_counts(A,B)
  iszero(iLa*iLb) && return _no_factor_slice_msg(A, B)

  Bc = make_concretetype(B{T})
  factors = [ :( $Bc($(Expr(:tuple,
                [:( indep_comp_getindex(p, $(la + iLa*(lb-1))) ) for lb in 1:iLb]...))) )
              for la in 1:iLa ]
  Expr(:tuple, factors...)
end

# The first factors of `p`, one per independent component of its second factor.
@generated function _first_factors(p::ProductTensorValue{A,B,T}) where {A,B,T}
  _, _, iLa, iLb = _prod_counts(A,B)
  iszero(iLa*iLb) && return _no_factor_slice_msg(A, B)

  Ac = make_concretetype(A{T})
  factors = [ :( $Ac($(Expr(:tuple,
                [:( indep_comp_getindex(p, $(la + iLa*(lb-1))) ) for la in 1:iLa]...))) )
              for lb in 1:iLb ]
  Expr(:tuple, factors...)
end

# The product tensor with first factor type `A` and second factors `bs`.
@generated function _product_from_second_factors(
  ::Type{A}, bs::NTuple{iLa,B}) where {A,iLa,B<:MultiValue}

  iLb = num_indep_components(B)
  comps = [ :( indep_comp_getindex(bs[$la], $lb) ) for lb in 1:iLb for la in 1:iLa ]
  :( ProductTensorValue{$A,$B}($(Expr(:tuple, comps...))) )
end

# A scalar valued action gives a `ScalarValue` second factor, so a product
# tensor maps to a product tensor whatever the action returns.
_product_from_second_factors(::Type{A}, bs::NTuple{iLa,B}) where {A,iLa,B<:Number} =
  ProductTensorValue{A,ScalarValue}(bs)

# The product tensor with second factor type `B` and first factors `as`.
@generated function _product_from_first_factors(
  ::Type{B}, as::NTuple{iLb,A}) where {B,iLb,A<:MultiValue}

  iLa = num_indep_components(A)
  comps = [ :( indep_comp_getindex(as[$lb], $la) ) for lb in 1:iLb for la in 1:iLa ]
  :( ProductTensorValue{$A,$B}($(Expr(:tuple, comps...))) )
end

"""
    _map_second_factor(f, p::ProductTensorValue)

The product tensor with the first factor of `p` and the images by `f` of its
second factors.

`p` is the sum of `aₗ ⊗ bₗ` over the independent components of its first factor,
and the result is the sum of the `aₗ ⊗ f(bₗ)`, that is the image of `p` by `f`
acting on its second index block alone.

`f` must be linear: a first factor with interdependent components only has the
second factors of its dependent components reproduced by linearity.

A scalar valued `f` gives a [`ScalarValue`](@ref) second factor, so the result
is a product tensor whatever `f` returns.
"""
# `@inline`: the `map` closure is emitted out of line otherwise, once per factor
# slice, which costs several times the whole mapping.
@inline _map_second_factor(f, p::ProductTensorValue{A}) where A =
  _product_from_second_factors(A, map(f, _second_factors(p)))

"""
    _map_first_factor(f, p::ProductTensorValue)

The product tensor with the images by `f` of the first factors of `p` and its
second factor.

This is the mirror of [`_map_second_factor`](@ref): `p` is the sum of
`aₗ ⊗ bₗ` over the independent components of its second factor, and the result
is the sum of the `f(aₗ) ⊗ bₗ`, the image of `p` by `f` acting on its first
index block alone. `f` must be linear for the same reason.
"""
# `@inline` for the same reason as [`_map_second_factor`](@ref).
@inline _map_first_factor(f, p::ProductTensorValue{A,B}) where {A,B} =
  _product_from_first_factors(B, map(f, _first_factors(p)))

