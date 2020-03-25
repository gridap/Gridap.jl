
#size(a<:MultiValue) = size(a.array)

#length(a<:MultiValue) = length(a.array)

function getindex(a::VectorValue,i::Integer)
    a.data[i]
end

function getindex(a::TensorValue,i::Integer)
    a.data[i]
end

function getindex(a::TensorValue{D1,D2},i::Integer,j::Integer) where {D1,D2}
    index = (j-1)*D1 + i
    a.data[index]
end

@inline iterate(a::MultiValue) = iterate(a.data)

@inline iterate(a::MultiValue, state) = iterate(a.data, state)

eachindex(a::MultiValue) = eachindex(a.data)

function CartesianIndices(a::MultiValue)
  CartesianIndices(get_array(a))
end

function LinearIndices(a::MultiValue)
  LinearIndices(get_array(a))
end
