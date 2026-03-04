struct VectorWithEntryRemoved{T,A} <: AbstractVector{T}
  a::A
  index::Int
  function VectorWithEntryRemoved(a::AbstractVector,index::Integer)
    A = typeof(a)
    T = eltype(a)
    @assert 1 <= index <= length(a)
    new{T,A}(a,index)
  end
end

Base.IndexStyle(::Type{<:VectorWithEntryRemoved}) = IndexLinear()

@propagate_inbounds function Base.getindex(v::VectorWithEntryRemoved, i::Integer)
  @boundscheck checkbounds(v,i)
  @inbounds i < v.index ? v.a[i] : v.a[i+1]
end

Base.size(v::VectorWithEntryRemoved) = (length(v.a)-1,)

function array_cache(v::VectorWithEntryRemoved)
  array_cache(v.a)
end

@propagate_inbounds function getindex!(cache,v::VectorWithEntryRemoved, i::Integer)
  @boundscheck checkbounds(v,i)
  @inbounds i < v.index ? getindex!(cache,v.a,i) : getindex!(cache,v.a,i+1)
end

function Base.sum(a::VectorWithEntryRemoved)
  sum(a.a) - a.a[a.index]
end

@propagate_inbounds function Base.setindex!(a::VectorWithEntryRemoved, v, i::Integer)
  @boundscheck checkbounds(a,i)
  @inbounds i < a.index ? (a.a[i] = v) : (a.a[i+1] = v)
end
