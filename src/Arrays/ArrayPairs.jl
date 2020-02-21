
"""
"""
function pair_arrays(a::AbstractArray,b::AbstractArray)
  ArrayPair(a,b)
end

function pair_arrays(a::CompressedArray,b::CompressedArray)
  if (a.ptrs === b.ptrs) || (a.ptrs == b.ptrs)
    values = [ (ai,bi) for (ai,bi) in zip(a.values,b.values) ]
    ptrs = a.ptrs
    return CompressedArray(values,ptrs)
  else
    return ArrayPair(a,b)
  end
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

"""
"""
function unpair_arrays(pair::AbstractArray{<:Tuple})
  p = testitem(pair)
  N = length(p)
  @assert N > 0
  _unpair_arrays(pair,Val{1}(),Val{N}())
end

function unpair_arrays(pair::ArrayPair)
  pair.a, pair.b
end

function  _unpair_arrays(pair,::Val{n},::Val{N}) where {n,N}
  pn = UnpairedArray(Val{n}(),pair)
  ps = _unpair_arrays(pair,Val{n+1}(),Val{N}())
  (pn, ps...)
end

function  _unpair_arrays(pair,::Val{N},::Val{N}) where N
  pn = UnpairedArray(Val{N}(),pair)
  (pn, )
end

struct UnpairedArray{I,T,N,P<:AbstractArray} <: AbstractArray{T,N}
  pair::P
  function UnpairedArray(::Val{I},pair::AbstractArray{<:Tuple,N}) where {I,N}
    p = testitem(pair)
    T = typeof(p[I])
    P = typeof(pair)
    new{I,T,N,P}(pair)
  end
end

Base.size(a::UnpairedArray) = size(a.pair)

Base.IndexStyle(::Type{UnpairedArray{I,T,N,P}}) where {T,N,I,P} = IndexStyle(P)

@inline function Base.getindex(a::UnpairedArray{I},i::Integer) where I
  a.pair[i][I]
end

@inline function Base.getindex(a::UnpairedArray{I,T,N},i::Vararg{Integer,N}) where {I,T,N}
  a.pair[i...][I]
end

uses_hash(::Type{<:UnpairedArray}) = Val{true}()

function array_cache(hash,p::UnpairedArray)
  array_cache(hash,p.pair)
end

@inline function getindex!(cache,a::UnpairedArray{I},i) where I
  p = getindex!(cache,a.pair,i)
  p[I]
end

