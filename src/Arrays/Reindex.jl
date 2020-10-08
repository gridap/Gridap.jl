"""
    reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
"""
struct Reindex{A} <: Map
  values::A
end

@inline return_type(k::Reindex,x...) = typeof(first(k.values))

@inline return_type(k::Reindex,x::AbstractArray...) = typeof(testitem(k,x...))

@inline return_cache(k::Reindex,f) = array_cache(k.values)

@inline function return_cache(k::Reindex,a::AbstractArray)
  gids = a
  vals = k.values
  T = eltype(vals)
  r = zeros(T,size(gids))
  c = CachedArray(r)
end

@inline evaluate!(cache,k::Reindex,i) = getindex!(cache,k.values,i)

@inline function evaluate!(cache,k::Reindex,gids::AbstractArray)
  c = cache
  setsize!(c,size(gids))
  r = c.array
  for i in eachindex(gids)
    @inbounds r[i] = k.values[gids[i]]
  end
  r
end

function reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
  lazy_map(Reindex(i_to_v),j_to_i)
end

function reindex(i_to_v::Fill, j_to_i::AbstractArray)
  v = i_to_v.value
  Fill(v,size(j_to_i)...)
end

function reindex(i_to_v::CompressedArray, j_to_i::AbstractArray)
  values = i_to_v.values
  ptrs = i_to_v.ptrs[j_to_i]
  CompressedArray(values,ptrs)
end

function reindex(gid_to_val::AbstractArray,lid_to_gid::AbstractArray{<:AbstractArray})
  lazy_map(Reindex(gid_to_val),lid_to_gid)
end

Base.getindex(a::LazyArray,i::AbstractArray) = reindex(a,i)

@propagate_inbounds function Base.setindex!(a::LazyArray{<:Fill{<:Reindex}},v,j::Integer)
  k = a.g.value
  i_to_v = k.values
  j_to_i, = a.f
  i = j_to_i[j]
  i_to_v[i]=v
end

@inline function testitem(k::Reindex,gids)
  gids == 0 ? testitem(k,1) :  evaluate(k,gids)#testvalue(eltype(k.values))
end
