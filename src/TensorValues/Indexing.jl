
getindex(arg::VectorValue,  i::Integer) = arg.data[i]
getindex(arg::VectorValue, ci::CartesianIndex{1}) = get_index(arg,ci[1])
getindex(arg::TensorValue,  i::Integer) = arg.data[i]

getindex(arg::TensorValue{D1,D2},i::Integer,j::Integer) where {D1,D2} = (index=((j-1)*D1)+i; arg.data[index])
getindex(arg::TensorValue{D1,D2},ci::CartesianIndex{2}) where {D1,D2} = get_index(arg,ci[1],ci[2])

@inline iterate(arg::MultiValue) = iterate(arg.data)
@inline iterate(arg::MultiValue, state) = iterate(arg.data, state)

eachindex(arg::MultiValue) = eachindex(arg.data)

CartesianIndices(arg::MultiValue) = CartesianIndices(get_array(arg))

LinearIndices(arg::MultiValue) = LinearIndices(get_array(arg))
