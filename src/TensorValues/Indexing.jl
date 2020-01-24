
size(a::MultiValue) = size(a.array)

length(a::MultiValue) = length(a.array)

@propagate_inbounds function getindex(
    a::MultiValue{S,T,N}, I::Vararg{Integer,N}) where {S,T,N}
    a.array[I...]
end

@propagate_inbounds function getindex(a::MultiValue, i::Integer)
    a.array[i]
end

eltype(a::Type{MultiValue{S,T,N,L}}) where {S,T,N,L} = T

@inline iterate(a::MultiValue) = iterate(a.array)

@inline iterate(a::MultiValue, state) = iterate(a.array, state)

eachindex(a::MultiValue) = eachindex(a.array)

function CartesianIndices(a::MultiValue)
  CartesianIndices(a.array)
end

function LinearIndices(a::MultiValue)
  LinearIndices(a.array)
end
