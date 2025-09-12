# It would be better to keep IndexLinear style for VectorValue, TensorValue and
# ThirdOrderTensorValue, their eachindex isn't as efficient as possible
Base.IndexStyle(::MultiValue) = IndexCartesian()
Base.IndexStyle(::Type{<:MultiValue}) = IndexCartesian()

# Necessary overloads due to wrong ::Number defaults
lastindex(arg::MultiValue) = length(arg)
lastindex(arg::MultiValue, d::Int) = (@inline; last(axes(arg, d)))

CartesianIndices(arg::MultiValue) = CartesianIndices(axes(arg))
LinearIndices(arg::MultiValue) = LinearIndices(axes(arg))

eachindex(arg::MultiValue) = CartesianIndices(arg)
eachindex(::IndexCartesian, arg::MultiValue) = eachindex(arg)
eachindex(::IndexLinear, arg::MultiValue) = SOneTo(length(arg))

keys(arg::MultiValue) = CartesianIndices(axes(arg))
keys(s::IndexStyle, arg::MultiValue) = eachindex(s, arg)

"""
    getindex(arg::MultiValue, inds...)
    getindex(arg::MultiValue, i::Integer)
    getindex(arg::MultiValue) = arg

`MultiValue`s can be indexed like `Base.AbstractArray`s.

The `Number` convention is used when no indices are provided: `arg[]` returns `arg`.

# Extented help

Similarly to StaticArray.jl, the type of the returned value depend on
inferability of the size from the type of the indices in `inds`.

Indeed, `getindex` returns
- a scalar (component) if `inds` contains only "scalars", i.e. `Integer` and `CartesianIndex`,
- a `<:MultiValue` tensor if `inds` additionally contains statically inferable index ranges/sets such as `Colon()`/`:`, `SOneTo` or `StaticArray`,
- an `Array` if `inds` additionally contains dynamically sized index sets such as `UnitRange`, `OneTo` or other types of `AbstractArray`.

# warning
  Indexing methods that do not return a scalars will loose any symmetry or other
  known internal constraint, as in `SymTensorValue(1:6...)[SOneTo(2),SOneTo(2)]`
""" # Those methods are necessary because MultiValue subtypes Number instead of AbstractArray
@propagate_inbounds getindex(arg::MultiValue, i::Integer) = getindex(arg, CartesianIndices(arg)[i])
@propagate_inbounds getindex(arg::MultiValue) = arg
# Size-inferable "scalar" indexing
const _ScalarIndices = Union{Integer, CartesianIndex}
@propagate_inbounds getindex(arg::MultiValue, inds::_ScalarIndices...) = getindex(arg, to_indices(arg, inds)...)
# Size-inferable "array" indexing,
const _StaticIndices = Union{Colon,SOneTo,StaticArray{<:Tuple, Int64},_ScalarIndices}
# the conversion to SArray is only necessary when there are CartesianIndex in `inds`, see https://github.com/JuliaArrays/StaticArrays.jl/issues/1059
@propagate_inbounds getindex(arg::MultiValue, inds::_StaticIndices...) = MultiValue(SArray(getindex(get_array(arg), inds...)))
# Not size-inferable "array" indexing, returns ::Base.Array. Includes ::OneTo, Array{Bool}, etc.
const _DynamicIndices= Union{AbstractArray,_StaticIndices}
@propagate_inbounds getindex(arg::MultiValue, inds::_DynamicIndices...) = getindex(get_array(arg), inds...)


# Cartesian indexing style implementation
@propagate_inbounds function getindex(arg::VectorValue,i::Integer)
  @boundscheck @check checkbounds(arg,i) === nothing
  @inbounds arg.data[i]
end

@propagate_inbounds function getindex(arg::TensorValue{D},i::Integer,j::Integer) where D
  @boundscheck @check checkbounds(arg,i,j) === nothing
  index = _2d_tensor_linear_index(D,i,j)
  @inbounds arg.data[index]
end

@propagate_inbounds function getindex(arg::AbstractSymTensorValue{D},i::Integer,j::Integer) where D
  @boundscheck @check checkbounds(arg,i,j) === nothing
  index = _2d_sym_tensor_linear_index(D,i,j)
  @inbounds arg.data[index]
end

@propagate_inbounds function getindex(arg::SymFourthOrderTensorValue{D},i::Integer,j::Integer,k::Integer,l::Integer) where D
  @boundscheck @check checkbounds(arg,i,j,k,l) === nothing
  index = _4d_sym_tensor_linear_index(D,i,j,k,l)
  @inbounds arg.data[index]
end

@propagate_inbounds function getindex(arg::ThirdOrderTensorValue{D1,D2},i::Integer,j::Integer,k::Integer) where {D1,D2}
  @boundscheck @check checkbounds(arg,i,j,k) === nothing
  index = _3d_tensor_linear_index(D1,D2,i,j,k)
  @inbounds arg.data[index]
end


function Base.checkbounds(A::MultiValue{S}, I::Integer...) where S
  if CartesianIndex(I...) ∉ CartesianIndices(A)
    throw(BoundsError(A,I))
  end
  nothing
end

@inline iterate(arg::MultiValue)        = iterate(arg.data)
@inline iterate(arg::MultiValue, state) = iterate(arg.data, state)

# This could be deprecated
data_index(::Type{<:VectorValue},i) = i
data_index(::Type{<:TensorValue{D}},i,j) where D = _2d_tensor_linear_index(D,i,j)
data_index(::Type{<:AbstractSymTensorValue{D}},i,j) where D = _2d_sym_tensor_linear_index(D,i,j)
data_index(::Type{<:ThirdOrderTensorValue{D1,D2}},i,j,k) where {D1,D2} = _3d_tensor_linear_index(D1,D2,i,j,k)
data_index(::Type{<:SymFourthOrderTensorValue{D}},i,j,k,l) where D = _4d_sym_tensor_linear_index(D,i,j,k,l)

_symmetric_index_gaps(i::Integer) = i*(i-1)÷2

_2d_tensor_linear_index(D,i,j) = ((j-1)*D)+i

_3d_tensor_linear_index(D1,D2,i,j,k) = (k-1)*D1*D2+(j-1)*D1+i

function _2d_sym_tensor_linear_index(D,i,j)
  _j,_i = minmax(i,j)
  index=_2d_tensor_linear_index(D,_i,_j)-_symmetric_index_gaps(_j)
  index
end

#function _4d_sym_tensor_linear_index(D,i,j,k,l)
#  _j=min(i,j)
#  _i=max(i,j)
#  _l=min(k,l)
#  _k=max(k,l)
#  block_length=_symmetric_index_gaps(D+1)
#  element_index=_2d_tensor_linear_index(D,_i,_j)-_symmetric_index_gaps(_j)
#  block_index=_2d_tensor_linear_index(D,_l,_k)-_symmetric_index_gaps(_l)
#  index=(block_index-1)*block_length+element_index
#  index
#end

function _4d_sym_tensor_linear_index(D,i,j,k,l)
  block_length = (D*(D+1))÷2
  block_index = _2d_sym_tensor_linear_index(D,i,j)
  element_index = _2d_sym_tensor_linear_index(D,k,l)
  index=(block_index-1)*block_length+element_index
  index
end
