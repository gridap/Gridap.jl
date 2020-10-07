# struct ArrayReindex{A} <: Mapping
  # values::A
# end

@inline return_type(k::Reindex,x::AbstractArray...) = typeof(testitem(k,x...))

@inline function return_cache(k::Reindex,a::AbstractArray)
  gids = a
  vals = k.values
  T = eltype(vals)
  r = zeros(T,size(gids))
  c = CachedArray(r)
end

@inline function evaluate!(cache,k::Reindex,gids::AbstractArray)
  c = cache
  setsize!(c,size(gids))
  r = c.array
  for i in eachindex(gids)
    @inbounds r[i] = k.values[gids[i]]
  end
  r
end

function reindex(gid_to_val::AbstractArray,lid_to_gid::AbstractArray{<:AbstractArray})
  lazy_map(Reindex(gid_to_val),lid_to_gid)
end

# function reindex(a::LazyArray{<:Fill{<:ArrayReindex}},b::AbstractArray)
#   map = a.g.value
#   gid_to_val = map.values
#   lid_to_gid, = a.f
#   lid_to_gid = reindex(lid_to_gid,b)
#   lazy_map(ArrayReindex(gid_to_val),lid_to_gid)
# end
