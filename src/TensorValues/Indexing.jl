
getindex(arg::MultiValue,  i::Integer) = arg.data[i]

function getindex(arg::TensorValue{D1,D2},i::Integer,j::Integer) where {D1,D2} 
    1 <= i <= D1 || throw(BoundsError(D1, i))
    1 <= j <= D2 || throw(BoundsError(D2, i))
    arg.data[_getindex(arg, i, j)]
end

function getindex(arg::SymTensorValue{D},i::Integer,j::Integer) where {D} 
    1 <= i <= D || throw(BoundsError(D, i))
    1 <= j <= D || throw(BoundsError(D, i))
    arg.data[_getindex(arg, i, j)]
end

getindex(arg::VectorValue, i::CartesianIndex{1}) = getindex(arg,ci[1])
getindex(arg::TensorValue{D1,D2},ci::CartesianIndex{2}) where {D1,D2} = getindex(arg,ci[1],ci[2])
getindex(arg::SymTensorValue{D1,D2},ci::CartesianIndex{2}) where {D1,D2} = getindex(arg,ci[1],ci[2])

@inline iterate(arg::MultiValue) = iterate(arg.data)
@inline iterate(arg::MultiValue, state) = iterate(arg.data, state)

eachindex(arg::MultiValue) = eachindex(arg.data)

CartesianIndices(arg::MultiValue) = CartesianIndices(get_array(arg))

LinearIndices(arg::MultiValue) = LinearIndices(get_array(arg))
