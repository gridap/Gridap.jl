# This Map has non-trivial domain, thus we need the define testargs
"""
    Reindex(values) -> Map
"""
struct Reindex{A} <: Map
  values::A
end

testargs(k::Reindex,i) = (one(i),)
testargs(k::Reindex,i::Integer...) = map(one,i)
return_cache(k::Reindex,i...) = array_cache(k.values)
@inline evaluate!(cache,k::Reindex,i...) = getindex!(cache,k.values,i...)

function testargs(k::Reindex,a::AbstractArray{T,N}) where {T,N}
  if length(a) == 0
    s = ntuple(i->1,Val{N}())
    b = similar(a,T,s)
    fill!(b,one(T))
    (b,)
  else
    (a,)
  end
end

function return_cache(k::Reindex,a::AbstractArray)
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

"""
    reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
"""
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

@propagate_inbounds function Base.setindex!(a::LazyArray{<:Fill{<:Reindex}},v,j::Integer)
  k = a.g.value
  i_to_v = k.values
  j_to_i, = a.f
  i = j_to_i[j]
  i_to_v[i]=v
end

