
eachindex(arg::MultiValue) = eachindex(1:prod(size(arg)))

CartesianIndices(arg::MultiValue) = CartesianIndices(size(arg))

LinearIndices(arg::MultiValue) = LinearIndices(size(arg))

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
