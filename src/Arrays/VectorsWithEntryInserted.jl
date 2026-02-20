mutable struct VectorWithEntryInserted{T,A} <: AbstractVector{T}
  a::A
  index::Int
  value::T
  function VectorWithEntryInserted(a::AbstractVector,index::Integer,value)
    A = typeof(a)
    T = eltype(a)
    @assert 1 <= index <= length(a)+1
    new{T,A}(a,index,value)
  end
end

Base.IndexStyle(::Type{<:VectorWithEntryInserted}) = IndexLinear()

@propagate_inbounds function Base.getindex(v::VectorWithEntryInserted, i::Integer)
  @boundscheck checkbounds(v,i)
  @inbounds i < v.index ? v.a[i] :
          ( i ==v.index ? v.value : v.a[i-1])
end

Base.size(v::VectorWithEntryInserted) = (length(v.a)+1,)

function array_cache(v::VectorWithEntryInserted)
  array_cache(v.a)
end

@propagate_inbounds function getindex!(cache,v::VectorWithEntryInserted, i::Integer)
  @boundscheck checkbounds(v,i)
  @inbounds i < v.index ? getindex!(cache,v.a,i) :
           (i ==v.index ? v.value : getindex!(cache,v.a,i-1))
end

function Base.sum(a::VectorWithEntryInserted)
  sum(a.a) + a.value
end

@propagate_inbounds function Base.setindex!(a::VectorWithEntryInserted, v, i::Integer)
  @boundscheck checkbounds(a,i)
  @inbounds i < a.index ? (a.a[i] = v) :
           (i ==a.index ? a.value = v : a.a[i-1]=v)
end
