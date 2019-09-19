module CachedSubVectors

using Base: @propagate_inbounds

export CachedSubVector
export locate!
import Base: size
import Base: getindex
import Base: setindex!
import Base: IndexStyle

mutable struct CachedSubVector{T,V<:AbstractArray{T,1}} <: AbstractArray{T,1}
  vector::V
  pini::Int
  pend::Int
end

function CachedSubVector(v::AbstractVector)
  CachedSubVector(v,1,length(v))
end

size(self::CachedSubVector) = (1+self.pend-self.pini,)

@propagate_inbounds function getindex(self::CachedSubVector, i::Int)
  @inbounds self.vector[self.pini+i-1]
end

@propagate_inbounds function setindex!(self::CachedSubVector,v, i::Int)
  @inbounds self.vector[self.pini+i-1] = v
end

IndexStyle(::Type{CachedSubVector{T,V}}) where {T,V} = IndexLinear()

function locate!(self::CachedSubVector,pini,pend)
  self.pini = pini
  self.pend = pend
end

end # module
