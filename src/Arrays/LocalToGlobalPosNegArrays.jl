
"""
"""
struct LocalToGlobalPosNegArray{T,M,N,L,V} <: AbstractArray{Array{T,M},N}
  lid_to_gid::L
  gid_to_val_pos::V
  gid_to_val_neg::V

  @doc """
  """
  function LocalToGlobalPosNegArray(
    lid_to_gid::AbstractArray{<:AbstractArray},
    gid_to_val_pos::AbstractArray,
    gid_to_val_neg::AbstractArray)

    L = typeof(lid_to_gid)
    V = typeof(gid_to_val_pos)
    T = eltype(V)
    M = ndims(eltype(L))
    N = ndims(L)
    new{T,M,N,L,V}(lid_to_gid,gid_to_val_pos,gid_to_val_neg)
  end
end

size(a::LocalToGlobalPosNegArray) = size(a.lid_to_gid)

function IndexStyle(::Type{LocalToGlobalPosNegArray{T,M,N,L,V}}) where {T,M,N,L,V}
  IndexStyle(L)
end

function array_cache(a::LocalToGlobalPosNegArray)
  gids = testitem(a.lid_to_gid)
  T = eltype(a.gid_to_val_pos)
  r = zeros(T,size(gids))
  c = CachedArray(r)
  cl = array_cache(a.lid_to_gid)
  (cl,c)
end

function getindex!(cache,a::LocalToGlobalPosNegArray,i::Integer...)
  (cl,c) = cache
  gids = getindex!(cl,a.lid_to_gid,i...)
  setsize!(c,size(gids))
  r = c.array
  for (i,gid) in enumerate(gids)
    if gid > 0
      r[i] = a.gid_to_val_pos[gid]
    elseif gid < 0
      r[i] = a.gid_to_val_neg[-gid]
    else
      @unreachable "Only positive or negative indices allowed, not zero."
    end
  end
  r
end

function getindex(a::LocalToGlobalPosNegArray,i::Integer...)
  cache = array_cache(a)
  getindex!(cache,a,i...)
end


