
eachindex(arg::MultiValue) = eachindex(1:prod(size(arg)))

lastindex(arg::MultiValue) = length(arg)

CartesianIndices(arg::MultiValue) = CartesianIndices(size(arg))

LinearIndices(arg::MultiValue) = LinearIndices(size(arg))

getindex(arg::VectorValue, i::Integer) = arg.data[i]

function getindex(arg::TensorValue{D},i::Integer,j::Integer) where D
  index = _2d_tensor_linear_index(D,i,j)
  arg.data[index]
end

function getindex(arg::SymTensorValue{D},i::Integer,j::Integer) where D
  index = _2d_sym_tensor_linear_index(D,i,j)
  arg.data[index]
end

function getindex(arg::SymFourthOrderTensorValue{D},i::Integer,j::Integer,k::Integer,l::Integer) where D
  index = _4d_sym_tensor_linear_index(D,i,j,k,l)
  arg.data[index]
end

function getindex(arg::ThirdOrderTensorValue{D1,D2},i::Integer,j::Integer,k::Integer) where {D1,D2}
  index = _3d_tensor_linear_index(D1,D2,i,j,k)
  arg.data[index]
end

getindex(arg::VectorValue, ci::CartesianIndex{1})              = getindex(arg,ci[1])
getindex(arg::TensorValue,ci::CartesianIndex{2})               = getindex(arg,ci[1],ci[2])
getindex(arg::SymTensorValue,ci::CartesianIndex{2})            = getindex(arg,ci[1],ci[2])
getindex(arg::ThirdOrderTensorValue,ci::CartesianIndex{3})     = getindex(arg,ci[1],ci[2],ci[3])
getindex(arg::SymFourthOrderTensorValue,ci::CartesianIndex{4}) = getindex(arg,ci[1],ci[2],ci[3],ci[4])

getindex(arg::MultiValue, i::Integer) = getindex(arg, CartesianIndices(arg)[i])
getindex(arg::TensorValue, i::Integer) = arg.data[i]
getindex(arg::ThirdOrderTensorValue, i::Integer) = arg.data[i]

@inline iterate(arg::MultiValue)        = iterate(arg.data)
@inline iterate(arg::MultiValue, state) = iterate(arg.data, state)

data_index(::Type{<:VectorValue},i) = i
data_index(::Type{<:TensorValue{D}},i,j) where D = _2d_tensor_linear_index(D,i,j)
data_index(::Type{<:SymTensorValue{D}},i,j) where D = _2d_sym_tensor_linear_index(D,i,j)
data_index(::Type{<:ThirdOrderTensorValue{D1,D2}},i,j,k) where {D1,D2} = _3d_tensor_linear_index(D1,D2,i,j,k)
data_index(::Type{<:SymFourthOrderTensorValue{D}},i,j,k,l) where D = _4d_sym_tensor_linear_index(D,i,j,k,l)

_symmetric_index_gaps(i::Integer) = i*(i-1)÷2

_2d_tensor_linear_index(D,i,j) = ((j-1)*D)+i

_3d_tensor_linear_index(D1,D2,i,j,k) = (k-1)*D1*D2+(j-1)*D1+i

function _2d_sym_tensor_linear_index(D,i,j)
  _j=min(i,j)
  _i=max(i,j)
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
