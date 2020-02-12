
"""
"""
function pair_arrays(a::AbstractArray,b::AbstractArray)
  ArrayPair(a,b)
end

struct ArrayPair{T,N,A,B} <:AbstractArray{T,N}
  a::A
  b::B

  function ArrayPair(a::AbstractArray,b::AbstractArray)
    @assert length(a) == length(b)
    A = typeof(a); B = typeof(b)
    Ta = eltype(a); Tb = eltype(b)
    T = Tuple{Ta,Tb}
    N = ndims(a)
    new{T,N,A,B}(a,b)
  end

end

Base.size(p::ArrayPair) = size(p.a)

Base.IndexStyle(::Type{ArrayPair{T,N,A,B}}) where {T,N,A,B} = IndexStyle(A)

@inline Base.getindex(p::ArrayPair,i::Integer) = (p.a[i], p.b[i])

@inline Base.getindex(p::ArrayPair{T,N},i::Vararg{Integer,N}) where {T,N} = (p.a[i...], p.b[i...])

uses_hash(::Type{<:ArrayPair}) = Val{true}()

function array_cache(hash,p::ArrayPair)
  ca = array_cache(hash,p.a)
  cb = array_cache(hash,p.b)
  (ca,cb)
end

@inline function getindex!(cache,p::ArrayPair,i)
  ca, cb = cache
  ai = getindex!(ca,p.a,i)
  bi = getindex!(cb,p.b,i)
  (ai,bi)
end

