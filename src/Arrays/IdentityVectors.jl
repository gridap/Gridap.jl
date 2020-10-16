
"""
"""
struct IdentityVector{T<:Integer} <: AbstractVector{T}
  length::T
end

function getindex(c::IdentityVector{T},i::Integer) where T
  @check i > 0
  @check i <= c.length
  j::T = i
  j
end

size(c::IdentityVector) = (c.length,)

IndexStyle(::Type{<:IdentityVector}) = IndexLinear()

function lazy_map(k::Reindex{<:AbstractArray}, indices::IdentityVector)
  @check length(k.values) == length(indices)
  k.values
end

#function lazy_map(k::Reindex{<:LazyArray},b::IdentityVector)
#  @check length(k.values) == length(indices)
#  k.values
#end

function lazy_map(k::Reindex{<:Fill},b::IdentityVector)
  @check length(k.values) == length(indices)
  k.values
end

function lazy_map(k::Reindex{<:CompressedArray},b::IdentityVector)
  @check length(k.values) == length(indices)
  k.values
end
