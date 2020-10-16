
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

