
"""
"""
struct BlockArrayCOO{T,N,A<:AbstractArray{T,N}} <: GridapType # TODO in the future this needs to implement block array
  blocks::Vector{A}
  coordinates::Vector{NTuple{N,Int}}
end
