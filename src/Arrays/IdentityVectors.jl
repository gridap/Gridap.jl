
"""
    identity_vector(l::Integer)
"""
function identity_vector(l::Integer)
  IdentityVector(l)
end

struct IdentityVector{T<:Integer} <: AbstractVector{T}
  length::T
end

function getindex(c::IdentityVector{T},i::Integer) where T
  @assert i > 0
  @assert i <= c.length
  j::T = i
  j
end

size(c::IdentityVector) = (c.length,)

IndexStyle(::Type{<:IdentityVector}) = IndexLinear()

function reindex(values::AbstractArray, indices::IdentityVector)
  @assert length(values) == length(indices)
  values
end

function reindex(a::AppliedArray,b::IdentityVector)
  a
end
