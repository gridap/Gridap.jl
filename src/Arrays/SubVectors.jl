
struct SubVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
  vector::A
  pini::Int
  pend::Int
end

Base.size(a::SubVector) = (1+a.pend-a.pini,)

@propagate_inbounds function Base.getindex(a::SubVector,i::Integer)
  a.vector[i+a.pini-1]
end

@propagate_inbounds function Base.setindex!(a::SubVector,v,i::Integer)
  a.vector[i+a.pini-1] = v
end

Base.IndexStyle(::Type{SubVector{T,A}}) where {T,A} = IndexLinear()

