###############################################################
# MultiValue Type
###############################################################

"""
    MultiValue{S,T,N,L} <: Number

Abstract type representing a multi-dimensional number value. The parameters are analog to that of StaticArrays.jl:
- `S` is a Tuple type holding the size of the tensor, e.g. Tuple{3} for a 3d vector or Tuple{2,4} for a 2 rows and 4 columns tensor,
- `T` is the type of the scalar components, should be subtype of `Number`,
- `N` is the order of the tensor, the length of `S`,
- `L` is the number of components stored internally.

`MultiValue`s are immutable.
"""
abstract type MultiValue{S,T,N,L} <: Number end

@inline Base.Tuple(arg::MultiValue) = arg.data

# Custom type printing

function show(io::IO,v::MultiValue)
  print(io,v.data)
end

function show(io::IO,::MIME"text/plain",v:: MultiValue)
  print(io,typeof(v))
  print(io,v.data)
end

###############################################################
# Introspection
###############################################################

eltype(::Type{<:MultiValue{S,T}}) where {S,T} = T
eltype(::MultiValue{S,T}) where {S,T} = T

length(::Type{<:MultiValue{S}}) where S = prod(size(MultiValue{S}))
length(a::MultiValue)  = length(typeof(a))

size(::Type{<:MultiValue{S}}) where S = _size_from_tuple_type(S)
size(a::MultiValue) = size(typeof(a))

_size_from_tuple_type(::Type{Tuple{}}) = ()
_size_from_tuple_type(::Type{Tuple{D}}) where D = (D,)
_size_from_tuple_type(::Type{Tuple{D1,D2}}) where {D1,D2} = (D1,D2)
_size_from_tuple_type(::Type{Tuple{D1,D2,D3}}) where {D1,D2,D3} = (D1,D2,D3)
_size_from_tuple_type(::Type{Tuple{D1,D2,D3,D4}}) where {D1,D2,D3,D4} = (D1,D2,D3,D4)
_size_from_tuple_type(::Type{S}) where S = @notimplemented

"""
    num_components(::Type{<:Number})
    num_components(a::Number)

Total number of components of a `Number` or `MultiValue`, that is 1 for scalars
and the product of the size dimensions for a `MultiValue`. This is the same as `length`.
"""
num_components(a::Number) = num_components(typeof(a))
num_components(::Type{<:Number}) = 1
num_components(T::Type{<:MultiValue}) = @unreachable "$T type is too abstract to count its components, the shape (firt parametric argument) is neede."
num_components(::Type{<:MultiValue{S}}) where S = length(MultiValue{S})

"""
    num_indep_components(::Type{<:Number})
    num_indep_components(a::Number)

Number of independent components of a `Number`, that is `num_components`
minus the number of components determined from others by symmetries or constraints.

For example, a `TensorValue{3,3}` has 9 independent components, a `SymTensorValue{3}`
has 6 and a `SymTracelessTensorValue{3}` has 5. But they all have 9 (non independent) components.
"""
num_indep_components(::Type{T}) where T<:Number = num_components(T)
num_indep_components(::T) where T<:Number = num_indep_components(T)

@deprecate n_components num_components


#######################################################################
# Other constructors and conversions implemented for more generic types
#######################################################################

zero(::Type{V}) where V<:MultiValue{S,T} where {S,T} = V(tfill(zero(T),Val(num_indep_components(V))))
zero(a::MultiValue) = zero(typeof(a))

one(::Type{V}) where V<:MultiValue = @unreachable "`one` not defined for $V"
one(a::MultiValue) = one(typeof(a))

function rand(rng::AbstractRNG,::Random.SamplerType{V}) where V<:MultiValue{D,T} where {D,T}
  Li = num_indep_components(V)
  vrand = rand(rng, SVector{Li,T})
  V(Tuple(vrand))
end

"""
    change_eltype(m::Number,::Type{T2})
    change_eltype(M::Type{<:Number},::Type{T2})

For multivalues, returns `M` or `typeof(m)` but with the component type (`MultiValue`'s parametric type `T`) changed to `T2`.

For scalars (or any non MultiValue number), `change_eltype` returns T2.
"""
change_eltype(::Type{<:Number}, T::Type) = T
change_eltype(::T,::Type{T2}) where {T<:Number,T2} = change_eltype(T,T2)

change_eltype(a::MultiValue,::Type{T2}) where T2 = change_eltype(typeof(a),T2)

get_array(a::MultiValue{S,T}) where {S,T} = convert(SArray{S,T},a)

"""
    Mutable(T::Type{<:MultiValue}) -> ::Type{<:MArray}
    Mutable(a::MultiValue)

Return the concrete mutable `MArray` type (defined by `StaticArrays.jl`) corresponding
to the `MultiValue` type T or array size and type of `a`.

See also [`mutable`](@ref).
"""
Mutable(::Type{<:MultiValue}) = @abstractmethod
# typeof(zero(...)) is for the N and L type parameters to be computed, and get
# e.g. MVector{D,T} == MMarray{Tuple{D},T,1,D} instead of MMarray{Tuple{D},T}
Mutable(::Type{<:MultiValue{S,T}}) where {S,T} = typeof(zero(MArray{S,T}))
Mutable(a::MultiValue) = Mutable(typeof(a))

"""
    mutable(a::MultiValue)

Converts `a` into a mutable array of type `MArray` defined by `StaticArrays.jl`.

See also [`Mutable`](@ref).
"""
mutable(a::MultiValue{S}) where S = MArray{S}(Tuple(get_array(a)))

###############################################################
# Conversions
###############################################################

# Direct conversion
convert(::Type{V}, arg::AbstractArray) where V<:MultiValue{S,T} where {S,T} = V(arg)
convert(::Type{V}, arg::Tuple) where V<:MultiValue{S,T} where {S,T} = V(arg)

# Inverse conversion
convert(::Type{<:NTuple{L,T}}, arg::MultiValue) where {L,T} = NTuple{L,T}(Tuple(arg))

###############################################################
# Indexing independant components
###############################################################

# This should probably not be exported, as (accessing) the data field of
# MultiValue is not a public api
"""
Transforms Cartesian indices to linear indices that index `MultiValue`'s private internal storage, this should'nt be used.
"""
function data_index(::Type{<:MultiValue},i...)
  @abstractmethod
end

# The order of export of components is that of their position in the .data
# field, but the actual method "choosing" the export order is
# Gridap.Visualization._prepare_data(::Multivalue).
"""
    indep_comp_getindex(a::Number,i)

Get the `i`th independent component of `a`. It only differs from `getindex(a,i)`
when the components of `a` are interdependent, see [`num_indep_components`](@ref).
`i` should be in `1:num_indep_components(a)`.
"""
@propagate_inbounds function indep_comp_getindex(a::Number,i)
  @boundscheck @check 1 <= i <= num_indep_components(Number)
  @inbounds a[i]
end

@propagate_inbounds function indep_comp_getindex(a::T,i) where {T<:MultiValue}
  @boundscheck @check 1 <= i <= num_indep_components(T)
  @inbounds _get_data(a,i)
end

# abstraction of Multivalue data access in case subtypes of MultiValue don't
# store its data in a data field
@propagate_inbounds function _get_data(a::MultiValue,i)
  a.data[i]
end

"""
    indep_components_names(::MultiValue)

Return an array of strings containing the component labels in the order they
are exported in VTK file.

If all dimensions of the tensor shape S are smaller than 3, the components
are named with letters "X","Y" and "Z" similarly to the automatic naming
of Paraview. Else, if max(S)>3, they are labeled by integers starting from "1".
"""
function indep_components_names(::Type{MultiValue{S,T,N,L}}) where {S,T,N,L}
  return ["$i" for i in 1:L]
end

@inline function ForwardDiff.value(x::MultiValue{S,<:ForwardDiff.Dual}) where S
  VT = change_eltype(x,ForwardDiff.valtype(eltype(x)))
  data = map(ForwardDiff.value,x.data)
  return VT(data)
end
