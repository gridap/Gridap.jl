###############################################################
# Abstract type for all symetric tensors
###############################################################

"""
    AbstractSymMultiValue{S,T,N,L} <: MultiValue{S,T,N,L}

Abstract type representing any fully symmetric tensor, that is a tensor whose
components are invariant under any permutation of its `N` indices:

    a[I...] == a[I[σ]...] for every permutation σ of 1:N

Equivalently, two components of indices `I` and `J` are equal if and only if
`sort(I) == sort(J)`, so a component is determined by the multiset of its indices.

See also [`AbstractSymTensorValue`](@ref), the second-order case.
"""
abstract type AbstractSymMultiValue{S,T,N,L} <: MultiValue{S,T,N,L} end


###############################################################
# SymNTensorValue Type
###############################################################

"""
    SymNTensorValue{S,T,N,L} <: AbstractSymMultiValue{S,T,N,L}

Type representing a fully symmetric `N`th-order tensor of size `S` and scalar
type `T`, storing its `L = binomial(D+N-1,N)` independent components, where `D`
is the dimension shared by all `N` indices.

Its components are invariant under any permutation of the indices, so `a[I...]
== a[J...]` if and only if `sort(I) == sort(J)`. Such a tensor holds for
instance the independent `N`th order derivatives of a function of `D`
variables.

It is constructed by providing the components of the non-decreasing indices
`1 ≤ i_1 ≤ … ≤ i_N ≤ D` in lexicographic order, e.g. for `N=3` and `D=2`:
`(1,1,1)`, `(1,1,2)`, `(1,2,2)`, `(2,2,2)`.

Use [`sorted_multisets`](@ref) to enumerate the indices of the independent
components, and [`multiset_index`](@ref) to locate one of them.

Each index carries its own variance, so `S` holds `Dim`s that are equal in
value but not necessarily in variance, e.g.

```julia
SymNTensorValue{Tuple{3,3,3}}(1:10...)         # un-varianced
SymNTensorValue{Tuple{con3,con3,con3}}(1:10...) # fully contravariant
SymNTensorValue{Tuple{cov2,con2}}(1,2,3)        # mixed variance
```
"""
struct SymNTensorValue{S,T,N,L} <: AbstractSymMultiValue{S,T,N,L}
  data::NTuple{L,T}

  function SymNTensorValue{S,T,N}(data::NTuple{L,T}) where {S,T,N,L}
    _, n, l = _sym_n_tensor_params(S)
    @assert n == N "wrong tensor order: got $N, expected $n"
    @assert l == L "wrong number of components: got $L, expected $l"
    new{S,T,N,L}(data)
  end
end

_sym_shape_error(::Type{S}) where S<:Tuple =
  "Invalid fully symmetric tensor size $S: it must hold `Dim`s that are all equal in value"

"""
    _sym_n_tensor_params(::Type{S}) -> (D, N, L)

Common dimension `D`, order `N` and number of independent components
`L = binomial(D+N-1,N)` of the fully symmetric tensor size `S`, erroring if `S`
is not a fully symmetric shape.
"""
function _sym_n_tensor_params(::Type{S}) where S<:Tuple
  p = _sym_shape(S)
  p === nothing && @unreachable _sym_shape_error(S)
  p
end

"""
    _sym_shape(::Type{S}) -> (D, N, L)

Validates `S` is a correct size Tuple for a fully symmetric tensor, and return
the common dimension `D` (un-varianced), it's order `N`, and number of
independent compenents `L`.
"""
function _sym_shape(::Type{S}) where S<:Tuple
  dims = size(MultiValue{S})
  all(d -> d isa Dim, dims) || return nothing
  ds = _rm_size_variance(dims)
  N = length(ds)
  N == 0 && return (0, 0, 1)
  D = first(ds)
  # `all(==(D), ds)` rather than `allequal(ds)`: only the former is provably
  # terminating, which is what lets this function be constant folded.
  all(==(D), ds) || return nothing
  (D, N, binomial(D+N-1, N))
end

function promote_rule(
  ::Type{<:SymNTensorValue{S,Ta,N}},
  ::Type{<:SymNTensorValue{S,Tb,N}}) where {S,Ta,Tb,N}

  T = promote_type(Ta,Tb)
  SymNTensorValue{S,T,N}
end

###############################################################
# Constructors (SymNTensorValue)
###############################################################

function SymNTensorValue{S,T}(data::NTuple{L,T}) where {S,T,L}
  _, N, _ = _sym_n_tensor_params(S)
  SymNTensorValue{S,T,N}(data)
end

# Empty constructors
SymNTensorValue{S}()                where S     = SymNTensorValue{S,Int}(NTuple{0,Int}())
SymNTensorValue{S,T}()              where {S,T} = SymNTensorValue{S,T}(NTuple{0,T}())
SymNTensorValue{S}(data::NTuple{0}) where S     = SymNTensorValue{S,Int}(data)

# NTuple argument constructors
SymNTensorValue{S}(data::NTuple{L,T}) where {S,L,T}               = SymNTensorValue{S,T}(data)
SymNTensorValue{S,T1}(data::NTuple{L,T2}) where {S,L,T1,T2}       = SymNTensorValue{S,T1}(NTuple{L,T1}(data))
SymNTensorValue{S,T1,N}(data::NTuple{L,T2}) where {S,L,T1,T2,N}   = SymNTensorValue{S,T1}(NTuple{L,T1}(data))
SymNTensorValue{S,T1,N,L}(data::NTuple{L,T2}) where {S,L,T1,T2,N} = SymNTensorValue{S,T1}(NTuple{L,T1}(data))

# single Tuple argument constructors
SymNTensorValue{S}(data::Tuple) where S               = SymNTensorValue{S}(promote(data...))
SymNTensorValue{S,T1}(data::Tuple) where {S,T1}       = SymNTensorValue{S,T1}(NTuple{length(data),T1}(data))
SymNTensorValue{S,T1,N}(data::Tuple) where {S,T1,N}   = SymNTensorValue{S,T1}(NTuple{length(data),T1}(data))
SymNTensorValue{S,T1,N,L}(data::Tuple) where {S,T1,N,L} = SymNTensorValue{S,T1}(NTuple{L,T1}(data))

# Vararg constructors
SymNTensorValue{S}(data::Number...) where S               = SymNTensorValue{S}(data)
SymNTensorValue{S,T1}(data::Number...) where {S,T1}       = SymNTensorValue{S,T1}(data)
SymNTensorValue{S,T1,N}(data::Number...) where {S,T1,N}   = SymNTensorValue{S,T1}(data)
SymNTensorValue{S,T1,N,L}(data::Number...) where {S,T1,N,L} = SymNTensorValue{S,T1}(data)

# single AbstractArray argument constructors

# Components of the non-decreasing indices, in the order they are stored. Like
# `_flatten_upper_triangle`, this reads `data` where it is symmetric and never
# checks that it is.
@generated function _flatten_sorted_components(data::AbstractArray, ::Type{S}) where S<:Tuple
  p = _sym_shape(S)
  p === nothing && return :( @unreachable $(_sym_shape_error(S)) )
  D, N, _ = p
  s = ntuple(_ -> D, N)
  msg = "expected an array of size $s"
  comps = [ :(data[$(I...)]) for I in sorted_multisets(D, N) ]
  quote
    Base.require_one_based_indexing(data)
    @check $s == size(data) $msg
    return $(Expr(:tuple, comps...))
  end
end

SymNTensorValue{S}(data::AbstractArray{T}) where {S,T}             = SymNTensorValue{S,T}(_flatten_sorted_components(data,S))
SymNTensorValue{S,T1}(data::AbstractArray{T2}) where {S,T1,T2}     = SymNTensorValue{S,T1}(_flatten_sorted_components(data,S))
SymNTensorValue{S,T1,N}(data::AbstractArray{T2}) where {S,T1,T2,N} = SymNTensorValue{S,T1}(_flatten_sorted_components(data,S))
SymNTensorValue{S,T1,N,L}(data::AbstractArray{T2}) where {S,T1,T2,N,L} = SymNTensorValue{S,T1}(_flatten_sorted_components(data,S))

###############################################################
# Conversions (SymNTensorValue)
###############################################################

# Full tensor components: `A[I] = a_{sort(I)}`.
#
# A `StaticArray` size parameter is a plain `Int`, so the variance of the
# dimensions is dropped here as it is by `get_array`.
@generated function _SymNTensorValue_to_array(arg::SymNTensorValue{Sa,T,N,L}) where {Sa,T,N,L}
  D, _, _ = _sym_shape(Sa) # validated by the constructor
  comps = Expr[]
  for idx in CartesianIndices(ntuple(_ -> D, N))
    l = multiset_index(sort!(collect(Tuple(idx))), D)
    push!(comps, :(arg.data[$l]))
  end

  S = NTuple{N,D}
  quote
    return SArray{$S,T,$N,$(D^N)}( tuple($(comps...)) )
  end
end

# Inverse conversion
convert(::Type{<:MArray{S,T}}, arg::SymNTensorValue) where {S,T} = MArray{S,T}(_SymNTensorValue_to_array(arg))
convert(::Type{<:SArray{S,T}}, arg::SymNTensorValue) where {S,T} = SArray{S,T}(_SymNTensorValue_to_array(arg))

# Internal conversion
convert(::Type{<:SymNTensorValue{S,T}}, arg::SymNTensorValue{S}) where {S,T} = SymNTensorValue{S,T}(Tuple(arg))
convert(::Type{<:SymNTensorValue{S,T}}, arg::SymNTensorValue{S,T}) where {S,T} = arg

###############################################################
# Other constructors and conversions (SymNTensorValue)
###############################################################

# Only orders 0 and 2 have a multiplicative neutral element: the scalar `one(T)`,
# and the Kronecker delta, which is fully symmetric and verifies `a⋅δ == a`.
@generated function one(::Type{<:SymNTensorValue{S,T}}) where {S,T}
  p = _sym_shape(S)
  p === nothing && return :( @unreachable $(_sym_shape_error(S)) )
  D, N, _ = p
  N == 0 && return :( one(T) )
  N == 2 || return :( @unreachable $("Fully symmetric tensors of order $N do not have a multiplicative neutral element.") )
  comps = [ i == j ? :(one(T)) : :(zero(T)) for (i,j) in sorted_multisets(D, N) ]
  :( return SymNTensorValue{S,T,2}( $(Expr(:tuple, comps...))) )
end

change_eltype(::Type{<:SymNTensorValue{S}},::Type{T2}) where {S,T2} = SymNTensorValue{S,T2}
change_eltype(::Type{SymNTensorValue{S,T1,N,L}},::Type{T2}) where {S,T1,T2,N,L} = SymNTensorValue{S,T2,N,L}

num_indep_components(::Type{<:SymNTensorValue{S}}) where S = _sym_n_tensor_params(S)[3]

###############################################################
# VTK export (SymNTensorValue)
###############################################################

function indep_components_names(::Type{<:SymNTensorValue{S}}) where S
  D, N, _ = _sym_n_tensor_params(S)
  if D > 3
    return [ join(I) for I in sorted_multisets(D,N) ]
  else
    c_name = ["X", "Y", "Z"]
    return [ join(c_name[i] for i in I) for I in sorted_multisets(D,N) ]
  end
end
