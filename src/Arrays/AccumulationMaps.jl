
# Note: 
# The following relies on Julia optimizing allocations when appending to 
# relatively small arrays. 
# It was created for usecases where the size of the output array varies a lot.
# If the size of the output array does not change much, it is probably better to use
# CachedArrays instead.

struct AppendEntriesMap{T,A} <: Map
  i_to_values :: A
  function AppendEntriesMap(i_to_values::AbstractArray{<:AbstractVector{T}}) where T
    A = typeof(i_to_values)
    new{T,A}(i_to_values)
  end
end

function return_cache(k::AppendEntriesMap{T}, I) where T
  c = array_cache(k.i_to_values)
  buffer = Vector{T}(undef,0) 
  return c, buffer
end

function evaluate!(cache, k::AppendEntriesMap{T}, I) where T
  c, buffer = cache
  empty!(buffer)
  for i in I
    v = getindex!(c, k.i_to_values, i)
    append!(buffer, v)
  end
  return buffer
end
