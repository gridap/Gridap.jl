
"""
"""
struct IdentityVector{T<:Integer} <: AbstractVector{T}
  length::T
end

@propagate_inbounds function getindex(c::IdentityVector{T},i::Integer) where T
  @boundscheck @check 0 < i <= c.length
  return T(i)
end

size(c::IdentityVector) = (c.length,)

IndexStyle(::Type{<:IdentityVector}) = IndexLinear()

