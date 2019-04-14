
mutable struct CachedArray{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  array::A
  size::NTuple{N,Int}
end

const CachedMatrix{T,A} = CachedArray{T,2,A}
const CachedVector{T,A} = CachedArray{T,1,A}

CachedArray(a::AbstractArray) = CachedArray(a,size(a))

size(self::CachedArray) = self.size

function setsize!(self::CachedArray{T,N},s::NTuple{N,Int}) where {T,N}
  @assert s <= size(self.array)
  self.size = s
end

@propagate_inbounds function getindex(self::CachedArray{T,N}, kj::Vararg{Integer,N}) where {T,N}
    @inbounds self.array[kj...]
end

@propagate_inbounds function setindex!(B::CachedArray{T,N}, v, kj::Vararg{Integer,N}) where {T,N}
    @inbounds B.array[kj...] = v
    v
end
