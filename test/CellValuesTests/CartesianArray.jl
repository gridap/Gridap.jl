
struct CartesianArray{T,N,A<:AbstractArray{T,N}} <:AbstractArray{T,N}
  data::A
end

Base.size(a::CartesianArray) = size(a.data)

function Base.getindex(
  a::CartesianArray{T,N},I::Vararg{Int,N}) where {T,N}
  a.data[I...]
end

