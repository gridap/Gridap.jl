
"""
    struct SkeletonPair{L,R} <: GridapType
      left::L
      right::R
    end
"""
struct SkeletonPair{L,R} <: GridapType
  left::L
  right::R
end

function reindex(a::AbstractArray,b::SkeletonPair)
  left = reindex(a,b.left)
  right = reindex(a,b.right)
  SkeletonPair(left,right)
end


