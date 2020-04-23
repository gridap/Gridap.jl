
eachindex(arg::MultiValue) = eachindex(1:prod(size(arg)))

CartesianIndices(arg::MultiValue) = CartesianIndices(size(arg))

LinearIndices(arg::MultiValue) = LinearIndices(size(arg))

_symmetric_index_gaps(i::Integer) = i*(i-1)÷2

_get_linear_index(D::Integer,i::Integer,j::Integer) = ((j-1)*D)+i

function _getindex(arg::TensorValue{D1,D2},i::Integer,j::Integer) where {D1,D2} 
    index=_get_linear_index(D1,i,j)
end

function _getindex(arg::SymTensorValue{D},i::Integer,j::Integer) where {D} 
    _j=min(i,j)
    _i=max(i,j)
    index=_get_linear_index(D,_i,_j)-_symmetric_index_gaps(_j)
end

function _getindex(arg::SymFourthOrderTensorValue{D},i::Integer,j::Integer,k::Integer,l::Integer) where {D}
    _j=min(i,j)
    _i=max(i,j)
    _l=min(k,l)
    _k=max(k,l)
    block_length=_symmetric_index_gaps(D+1)
    element_index=_get_linear_index(D,_i,_j)-_symmetric_index_gaps(_j)
    block_index=_get_linear_index(D,_l,_k)-_symmetric_index_gaps(_l)
    index=(block_index-1)*block_length+element_index
end

getindex(arg::VectorValue, i::Integer)                                                = arg.data[i]
getindex(arg::TensorValue,i::Integer,j::Integer)                                      = arg.data[_getindex(arg, i, j)]
getindex(arg::SymTensorValue,i::Integer,j::Integer)                                   = arg.data[_getindex(arg, i, j)]
getindex(arg::SymFourthOrderTensorValue,i::Integer,j::Integer,k::Integer,l::Integer)  = arg.data[_getindex(arg, i, j, k, l)]

getindex(arg::VectorValue, ci::CartesianIndex{1})              = getindex(arg,ci[1])
getindex(arg::TensorValue,ci::CartesianIndex{2})               = getindex(arg,ci[1],ci[2])
getindex(arg::SymTensorValue,ci::CartesianIndex{2})            = getindex(arg,ci[1],ci[2])
getindex(arg::SymFourthOrderTensorValue,ci::CartesianIndex{4}) = getindex(arg,ci[1],ci[2],ci[3],ci[4])

getindex(arg::MultiValue, i::Integer) = getindex(arg, CartesianIndices(arg)[i])

@inline iterate(arg::MultiValue)        = iterate(arg.data)
@inline iterate(arg::MultiValue, state) = iterate(arg.data, state)
