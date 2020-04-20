
getindex(arg::MultiValue, i::Integer) = arg.data[i]

function getindex(arg::TensorValue{D1,D2},i::Integer,j::Integer) where {D1,D2} 
    1 <= i <= D1 || throw(BoundsError(D1, i))
    1 <= j <= D2 || throw(BoundsError(D2, j))
    getindex(arg, _getindex(arg, i, j))
end

function getindex(arg::SymTensorValue{D},i::Integer,j::Integer) where {D} 
    1 <= i <= D || throw(BoundsError(D, i))
    1 <= j <= D || throw(BoundsError(D, j))
    getindex(arg, _getindex(arg, i, j))
end

function getindex(arg::SymFourthOrderTensorValue{D},i::Integer,j::Integer,k::Integer,l::Integer) where {D} 
    1 <= i <= D || throw(BoundsError(D, i))
    1 <= j <= D || throw(BoundsError(D, j))
    1 <= k <= D || throw(BoundsError(D, k))
    1 <= l <= D || throw(BoundsError(D, l))
    getindex(arg, _getindex(arg, i, j, k, l))
end

getindex(arg::VectorValue, ci::CartesianIndex{1}) = getindex(arg,ci[1])
getindex(arg::TensorValue{D1,D2},ci::CartesianIndex{2}) where {D1,D2} = getindex(arg,ci[1],ci[2])
getindex(arg::SymTensorValue{D},ci::CartesianIndex{2}) where {D} = getindex(arg,ci[1],ci[2])
getindex(arg::SymFourthOrderTensorValue{D},ci::CartesianIndex{4}) where {D} = getindex(arg,ci[1],ci[2],ci[3],ci[4])

@inline iterate(arg::MultiValue) = iterate(arg.data)
@inline iterate(arg::MultiValue, state) = iterate(arg.data, state)

eachindex(arg::MultiValue) = eachindex(arg.data)

CartesianIndices(arg::MultiValue) = CartesianIndices(size(arg))

LinearIndices(arg::MultiValue) = LinearIndices(size(arg))
