"""
    lazy_map(f,a::AbstractArray...) -> AbstractArray

Applies the `Map` (or `Function`) `f` to the entries of the arrays in `a`
(see the definition of [`Map`](@ref)).

The resulting array `r` is such that `r[i]` equals to `evaluate(f,ai...)` where `ai`
is the tuple containing the `i`-th entry of the arrays in `a` (see function
[`evaluate`](@ref) for more details).
In other words, the resulting array is numerically equivalent to:

    map( (x...)->evaluate(f,x...), a...)

# Examples

Using a function as mapping

```jldoctest
using Gridap.Arrays

a = collect(0:5)
b = collect(10:15)

c = lazy_map(+,a,b)

println(c)

# output
[10, 12, 14, 16, 18, 20]
```

Using a user-defined mapping

```jldoctest
using Gridap.Arrays
import Gridap.Arrays: evaluate!

a = collect(0:5)
b = collect(10:15)

struct MySum <: Map end

evaluate!(cache,::MySum,x,y) = x + y

k = MySum()

c = lazy_map(k,a,b)

println(c)

# output
[10, 12, 14, 16, 18, 20]
```
"""
function lazy_map(k,f::AbstractArray...)
  fi = map(testitem,f)
  T = return_type(k, fi...)
  lazy_map(k,T,f...)
end

#@inline lazy_map(::typeof(evaluate),k::AbstractArray,f::AbstractArray...) = LazyArray(k,f...)

# This is the function to be overload to specialize on the Map f
"""
    lazy_map(f,::Type{T},a::AbstractArray...) where T

Like [`lazy_map(f,a::AbstractArray...)`](@ref), but the user provides the element type
of the resulting array in order to circumvent type inference.
"""
@inline function lazy_map(k,T::Type,f::AbstractArray...)
  s = _common_size(f...)
  lazy_map(evaluate,T,Fill(k, s), f...)
end

# This is the function to be overload to specialize on the array types
@inline function lazy_map(::typeof(evaluate),T::Type,k::AbstractArray,f::AbstractArray...)
  s = _common_size(k,f...)
  N = length(s)
  LazyArray(T,Val(N),k,f...)
end

"""
Subtype of `AbstractArray` which is the result of `lazy_map`. It represents the
result of lazy_maping a `Map` to a set of arrays that
contain the mapping arguments. This struct makes use of the cache provided
by the mapping in order to compute its indices (thus allowing to prevent
allocation). The array is lazy, i.e., the values are only computed on
demand. It extends the `AbstractArray` API with two methods:

   `array_cache(a::AbstractArray)`
   `getindex!(a::AbstractArray,i...)`
"""
struct LazyArray{G,T,N,F} <: AbstractArray{T,N}
  g::G
  f::F
  function LazyArray(::Type{T}, g::AbstractArray, f::AbstractArray...) where T
    G = typeof(g)
    F = typeof(f)
    N = ndims(g)
    new{G,T,N,F}(g, f)
  end
  function LazyArray(::Type{T},::Val{N}, g::AbstractArray, f::AbstractArray...) where {T,N}
    @check ndims(g) == N || N == 1
    G = typeof(g)
    F = typeof(f)
    new{G,T,N,F}(g, f)
  end
end

#function LazyArray(g::AbstractArray{S}, f::AbstractArray...) where S
#  isconcretetype(S) ? gi = testitem(g) : @notimplemented
#  fi = map(testitem,f)
#  T = return_type(gi, fi...)
#  LazyArray(T, g, f...)
#end

IndexStyle(::Type{<:LazyArray}) = IndexCartesian()

IndexStyle(::Type{<:LazyArray{G,T,1} where {G,T}}) = IndexLinear()

uses_hash(::Type{<:LazyArray}) = Val{true}()

function array_cache(hash::Dict,a::LazyArray)
  function _getid(hash,cache::T,id) where T
    value::T = hash[id]
    value
  end
  id = objectid(a)
  cache = _array_cache!(hash,a)
  if ! haskey(hash,id)
    hash[id] = cache
  end
  _getid(hash,cache,id)
end

mutable struct IndexItemPair{T,V}
  index::T
  item::V
end

function _array_cache!(hash::Dict,a::LazyArray)
  @boundscheck begin
    @notimplementedif ! all(map(isconcretetype, map(eltype, a.f)))
    if ! (eltype(a.g) <: Function)
      @notimplementedif ! isconcretetype(eltype(a.g))
    end
  end
  gi = testitem(a.g)
  fi = map(testitem,a.f)
  cg = array_cache(hash,a.g)
  cf = map(fi->array_cache(hash,fi),a.f)
  cgi = return_cache(gi, fi...)
  index = -1
  #item = evaluate!(cgi,gi,testargs(gi,fi...)...)
  item = return_value(gi,fi...)
  (cg, cgi, cf), IndexItemPair(index, item)
end

@inline getindex!(cache, a::LazyArray, i::Integer) = _getindex_optimized!(cache,a,i)

@inline function getindex!(cache, a::LazyArray{G,T,N}, i::Vararg{Integer,N}) where {G,T,N}
  _getindex_optimized!(cache,a,i...)
end

@inline function _getindex_optimized!(cache, a::LazyArray, i...)
  _cache, index_and_item = cache
  index = LinearIndices(a)[i...]
  if index_and_item.index != index
    item = _getindex!(_cache,a,i...)
    index_and_item.index = index
    index_and_item.item = item
  end
  index_and_item.item
end

@inline function _getindex!(cache, a::LazyArray, i...)
  cg, cgi, cf = cache
  gi = getindex!(cg, a.g, i...)
  fi = map((cj,fj) -> getindex!(cj,fj,i...),cf,a.f)
  vi = evaluate!(cgi, gi, fi...)
  vi
end

function Base.getindex(a::LazyArray, i::Integer)
  gi = a.g[i]
  fi = map(fj -> fj[i],a.f)
  vi = evaluate(gi, fi...)
  vi
end

function Base.getindex(a::LazyArray{G,T,N}, i::Vararg{Integer,N}) where {G,T,N}
  gi = a.g[i...]
  fi = map(fj -> fj[i...],a.f)
  vi = evaluate(gi, fi...)
  vi
end

Base.size(a::LazyArray) = size(a.g)
Base.size(a::LazyArray{G,T,1} where {G,T}) = (length(a.g),)

function Base.sum(a::LazyArray)
  cache = array_cache(a)
  _sum_lazy_array(cache,a)
end

function _sum_lazy_array(cache,a)
  r = zero(eltype(a))
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
    r += ai
  end
  r
end

function testitem(a::LazyArray{A,T} where A) where T
  if length(a) > 0
    first(a)
  else
    gi = testitem(a.g)
    fi = map(testitem,a.f)
    return_value(gi,fi...)
  end::T
end

# Particular implementations for Fill

#function lazy_map(::typeof(evaluate),f::Fill, a::Fill...)
#  ai = map(ai->ai.value,a)
#  r = evaluate(f.value, ai...)
#  s = _common_size(f, a...)
#  Fill(r, s)
#end

function lazy_map(::typeof(evaluate),::Type{T}, f::Fill, a::Fill...) where T
  #lazy_map(evaluate, f, a...)
  ai = map(ai->ai.value,a)
  r = evaluate(f.value, ai...)
  s = _common_size(f, a...)
  Fill(r, s)
end

function _common_size(a::AbstractArray...)
  a1, = a
  @check all(map(ai->length(a1) == length(ai),a)) "Array sizes $(map(size,a)) are not compatible."
  if all( map(ai->size(a1) == size(ai),a) )
    size(a1)
  else
    (length(a1),)
  end
end

# Needed only for testing purposes
struct ArrayWithCounter{T,N,A} <: AbstractArray{T,N}
  array::A
  counter::Array{Int,N}
  function ArrayWithCounter(a::AbstractArray{T,N}) where {T,N}
    c = zeros(Int,size(a))
    new{T,N,typeof(a)}(a,c)
  end
end

Base.size(a::ArrayWithCounter) = size(a.array)

function Base.getindex(a::ArrayWithCounter,i::Integer)
  a.counter[i] += 1
  a.array[i]
end

function Base.getindex(a::ArrayWithCounter{T,N},i::Vararg{Integer,N}) where {T,N}
  a.counter[i...] += 1
  a.array[i...]
end

Base.IndexStyle(::Type{<:ArrayWithCounter{T,N,A}}) where {T,A,N} = IndexStyle(A)

function resetcounter!(a::ArrayWithCounter)
  fill!(a.counter,zero(eltype(a.counter)))
end

