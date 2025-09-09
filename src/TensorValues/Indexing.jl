
eachindex(arg::MultiValue) = eachindex(1:prod(size(arg)))
lastindex(arg::MultiValue) = length(arg)

CartesianIndices(arg::MultiValue) = CartesianIndices(size(arg))
LinearIndices(arg::MultiValue) = LinearIndices(size(arg))

Base.IndexStyle(::MultiValue) = IndexCartesian()
Base.IndexStyle(::Type{<:MultiValue}) = IndexCartesian()
Base.keys(::IndexLinear, A::MultiValue) = Base.OneTo(length(A))

"""
    getindex(arg::MultiValue, inds...)
    getindex(arg::MultiValue, i::Integer)
    getindex(arg::MultiValue) = arg

`MultiValue`s can be indexed like `Base.Array`s.

When indexing using a tupple `inds` (`IndexCartesian` style), `inds` may contain
`Integer`s and `CartesianIndex`s , but no `S/MArray`s (unlike arrays of `StaticArrays`).

The `Number` convention is used when no indices are provided: `arg` is returned.
""" # Those three methods are necessary because MultiValue subtypes Number instead of AbstractArray
@propagate_inbounds getindex(arg::MultiValue, inds...) = getindex(arg, to_indices(arg, inds)...)
@propagate_inbounds getindex(arg::MultiValue, i::Integer) = getindex(arg, CartesianIndices(arg)[i])
@propagate_inbounds getindex(arg::MultiValue) = @propagate_inbounds arg

# Slice getindex (1D)
@propagate_inbounds function getindex(arg::MultiValue, r::UnitRange{Int})
  @boundscheck @check checkbounds(arg,r) === nothing
  ntuple(i -> arg.data[r.start + (i-1)], length(r))
end

@propagate_inbounds function Base.getindex(arg::MultiValue, ::Colon)
  ntuple(i -> arg.data[i], length(arg))
end

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

function Base.checkbounds(A::MultiValue{S}, r::UnitRange{Int}) where S
  if any(i -> i < 1 || i > length(A), r)
    throw(BoundsError(A,r))
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
