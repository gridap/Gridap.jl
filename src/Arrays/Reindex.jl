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

"""
    reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
"""
function reindex(i_to_v::AbstractArray, j_to_i::AbstractArray)
  lazy_map(Reindex(i_to_v),j_to_i)
end

function lazy_map(k::Reindex{<:Fill}, j_to_i::AbstractArray)
  v = k.values.value
  Fill(v,size(j_to_i)...)
end

function lazy_map(k::Reindex{<:CompressedArray}, j_to_i::AbstractArray)
  i_to_v = k.values
  values = i_to_v.values
  ptrs = lazy_map(Reindex(i_to_v.ptrs),j_to_i)
  CompressedArray(values,ptrs)
end

@propagate_inbounds function Base.setindex!(a::LazyArray{<:Fill{<:Reindex}},v,j::Integer)
  k = a.g.value
  i_to_v = k.values
  j_to_i, = a.f
  i = j_to_i[j]
  i_to_v[i]=v
end

