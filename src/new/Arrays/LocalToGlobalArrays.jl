

"""
"""
struct LocalToGlobalArray{T,M,N,L,V} <: AbstractArray{Array{T,M},N}
  lid_to_gid::L
  gid_to_val::V

  @doc """
  """
  function LocalToGlobalArray(
    lid_to_gid::AbstractArray{<:AbstractArray},
    gid_to_val::AbstractArray)

    L = typeof(lid_to_gid)
    V = typeof(gid_to_val)
    T = eltype(V)
    M = ndims(eltype(L))
    N = ndims(L)
    new{T,M,N,L,V}(lid_to_gid,gid_to_val)
  end
end

size(a::LocalToGlobalArray) = size(a.lid_to_gid)

function IndexStyle(::Type{LocalToGlobalArray{T,M,N,L,V}}) where {T,M,N,L,V}
  IndexStyle(L)
end

function array_cache(a::LocalToGlobalArray)
  gids = testitem(a.lid_to_gid)
  T = eltype(a.gid_to_val)
  r = zeros(T,size(gids))
  c = CachedArray(r)
  cl = array_cache(a.lid_to_gid)
  (cl,c)
end

function getindex!(cache,a::LocalToGlobalArray,i::Integer...)
  (cl,c) = cache
  gids = getindex!(cl,a.lid_to_gid,i...)
  setsize!(c,size(gids))
  r = c.array
  for i in eachindex(gids)
    @inbounds r[i] = a.gid_to_val[gids[i]]
  end
  r
end

function getindex(a::LocalToGlobalArray,i::Integer...)
  cache = array_cache(a)
  getindex!(cache,a,i...)
end


