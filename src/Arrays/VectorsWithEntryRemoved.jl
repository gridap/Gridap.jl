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

function Base.getindex(v::VectorWithEntryRemoved,i::Integer)
  i < v.index ? v.a[i] : v.a[i+1]
end

Base.size(v::VectorWithEntryRemoved) = (length(v.a)-1,)

function array_cache(v::VectorWithEntryRemoved)
  array_cache(v.a)
end

function getindex!(cache,v::VectorWithEntryRemoved,i::Integer)
  i < v.index ? getindex!(cache,v.a,i) : getindex!(cache,v.a,i+1)
end

function Base.sum(a::VectorWithEntryRemoved)
  sum(a.a) - a.a[a.index]
end

function Base.setindex!(a::VectorWithEntryRemoved,v,i::Integer)
  i < a.index ? (a.a[i] = v) : (a.a[i+1] = v)
end
