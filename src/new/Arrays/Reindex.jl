
"""
    reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
"""
function reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
  Reindexed(i_to_v,j_to_i)
end

function reindex(i_to_v::Fill, j_to_i::AbstractArray)
  v = i_to_v.value
  l = length(j_to_i)
  Fill(v,l)
end

function reindex(i_to_v::CompressedArray, j_to_i::AbstractArray)
  values = i_to_v.values
  ptrs = i_to_v.ptrs[j_to_i]
  CompressedArray(values,ptrs)
end

struct Reindexed{T,N,A,B} <: AbstractArray{T,N}
  i_to_v::A
  j_to_i::B
  function Reindexed(i_to_v::AbstractArray, j_to_i::AbstractArray)
    T = eltype(i_to_v)
    N = ndims(j_to_i)
    A = typeof(i_to_v)
    B = typeof(j_to_i)
    new{T,N,A,B}(i_to_v,j_to_i)
  end
end

Base.size(a::Reindexed) = size(a.j_to_i)

Base.IndexStyle(::Type{Reindexed{T,N,A,B}}) where {T,N,A,B} = IndexStyle(B)

@propagate_inbounds function getindex(a::Reindexed,j::Integer)
  i = a.j_to_i[j]
  a.i_to_v[i]
end

@propagate_inbounds function getindex(a::Reindexed{T,N}, j::Vararg{Int,N}) where {T,N}
  i = a.j_to_i[j...]
  a.i_to_v[i]
end

array_cache(a::Reindexed) = array_cache(a.i_to_v)

function getindex!(cache,a::Reindexed,j::Integer...)
  i = a.j_to_i[j...]
  getindex!(cache,a.i_to_v,i)
end
