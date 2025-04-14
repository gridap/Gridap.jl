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

#lazy_map(::typeof(evaluate),k::AbstractArray,f::AbstractArray...) = LazyArray(k,f...)

# This is the function to be overload to specialize on the Map f
"""
    lazy_map(f,::Type{T},a::AbstractArray...) where T

Like [`lazy_map(f,a::AbstractArray...)`](@ref), but the user provides the element type
of the resulting array in order to circumvent type inference.
"""
function lazy_map(k,T::Type,f::AbstractArray...)
  s = _common_size(f...)
  lazy_map(evaluate,T,Fill(k, s), f...)
end

# This is the function to be overload to specialize on the array types
function lazy_map(::typeof(evaluate),T::Type,k::AbstractArray,f::AbstractArray...)
  s = _common_size(k,f...)
  N = length(s)
  LazyArray(T,Val(N),k,f...)
end

"""
Subtype of `AbstractArray` which is the result of `lazy_map`. It represents the
result of lazy_mapping a `Map` to a set of arrays that
contain the mapping arguments. This struct makes use of the cache provided
by the mapping in order to compute its indices (thus allowing to prevent
allocation). The array is lazy, i.e., the values are only computed on
demand. It extends the `AbstractArray` API with two methods:

   `array_cache(a::AbstractArray)`
   `getindex!(cache,a::AbstractArray,i...)`
"""
struct LazyArray{G,T,N,F} <: AbstractArray{T,N}
  maps::G
  args::F
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

function same_branch(a,b)
  a === b
end

function same_branch(a::Fill,b::Fill)
  typeof(a) != typeof(b) && return false
  size(a) != size(b) && return false
  a.value == b.value
end

function all_same_branch(a::Tuple,b::Tuple)
  for i in 1:length(a)
    if same_branch(a[i],b[i]) == false
      return false
    end
  end
  true
end

function same_branch(a::LazyArray,b::LazyArray)
  typeof(a) != typeof(b) && return false
  length(a.args) != length(b.args) && return false
  same_branch(a.maps,b.maps) && all_same_branch(a.args,b.args)
end

function _get_cache(dict,a)
  id = objectid(a)
  if haskey(dict,id)
    o,c = dict[id]
    return c
  end
  for item in dict
    if same_branch(a,second(item)[1])
      return second(item)[2]
    end
  end
  return nothing
end

function array_cache(dict::Dict,a::LazyArray)
  cache = _get_cache(dict,a)
  if cache === nothing
    _cache = _array_cache!(dict,a)
    dict[objectid(a)] = (a,_cache)
  else
    _cache = cache
  end
  _cache
end

mutable struct IndexItemPair{T,V}
  index::T
  item::V
end

function _array_cache!(dict::Dict,a::LazyArray)
  @boundscheck begin
    if ! all(map(isconcretetype, map(eltype, a.args)))
      for n in 1:length(a.args)
        @notimplementedif ! all(map(isconcretetype, map(eltype, a.args[n])))
      end
    end
    if ! (eltype(a.maps) <: Function)
      @notimplementedif ! isconcretetype(eltype(a.maps))
    end
  end
  gi = testitem(a.maps)
  fi = map(testitem,a.args)
  cg = array_cache(dict,a.maps)
  cf = map(fi->array_cache(dict,fi),a.args)
  cgi = return_cache(gi, fi...)
  index = -1
  #item = evaluate!(cgi,gi,testargs(gi,fi...)...)
  item = return_value(gi,fi...)
  (cg, cgi, cf), IndexItemPair(index, item)
end

function getindex!(cache, a::LazyArray, i::Integer)
  _cache, index_and_item = cache
  index = LinearIndices(a)[i]
  if index_and_item.index != index
    cg, cgi, cf = _cache
    gi = getindex!(cg, a.maps, i)
    index_and_item.item = _getindex_and_call!(cgi,gi,cf,a.args,i)
    index_and_item.index = index
  end
  index_and_item.item
end

function getindex!(cache, a::LazyArray{G,T,N}, i::Vararg{Integer,N}) where {G,T,N}
  _cache, index_and_item = cache
  index = LinearIndices(a)[i...]
  if index_and_item.index != index
    cg, cgi, cf = _cache
    gi = getindex!(cg, a.maps, i...)
    index_and_item.item = _getindex_and_call!(cgi,gi,cf,a.args,i...)
    index_and_item.index = index
  end
  index_and_item.item
end

function _getindex_and_call!(cgi,gi,cf,args,i...)
  fi = map((cj,fj) -> getindex!(cj,fj,i...),cf,args)
  evaluate!(cgi, gi, fi...)
end

function Base.getindex(a::LazyArray, i::Integer)
  gi = a.maps[i]
  fi = map(fj -> fj[i],a.args)
  vi = evaluate(gi, fi...)
  vi
end

function Base.getindex(a::LazyArray{G,T,N}, i::Vararg{Integer,N}) where {G,T,N}
  gi = a.maps[i...]
  fi = map(fj -> fj[i...],a.args)
  vi = evaluate(gi, fi...)
  vi
end

Base.size(a::LazyArray) = size(a.maps)
Base.size(a::LazyArray{G,T,1} where {G,T}) = (length(a.maps),)

function Base.sum(a::LazyArray)
  cache = array_cache(a)
  _sum_lazy_array(cache,a)
end

function _sum_lazy_array(cache,a)
  r = zero(testitem(a))
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
    gi = testitem(a.maps)
    fi = map(testitem,a.args)
    return_value(gi,fi...)
  end::T
end

# Particular implementations for Fill

function lazy_map(::typeof(evaluate),f::Fill, a::Fill...)
  ai = map(ai->ai.value,a)
  r = evaluate(f.value, ai...)
  s = _common_size(f, a...)
  Fill(r, s)
end

function lazy_map(::typeof(evaluate),::Type{T}, f::Fill, a::Fill...) where T
  ai = map(ai->ai.value,a)
  r = evaluate(f.value, ai...)
  s = _common_size(f, a...)
  Fill(r, s)
end

function _common_size(a::AbstractArray...)
  a1, = a
  #@check all(map(ai->length(a1) == length(ai),a)) "Array sizes $(map(size,a)) are not compatible."
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


# These extra methods are introduced to circumvent an unwanted run-time dispatch with julia 0.15 and 0.16
# see https://discourse.julialang.org/t/performance-depends-dramatically-on-compilation-order/58425
# Hopefully, they can be removed in the future

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  evaluate!(cgi,gi,f1)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  evaluate!(cgi,gi,f1,f2)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  f3 = getindex!(cf[3],args[3],i...)
  evaluate!(cgi,gi,f1,f2,f3)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any,Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  f3 = getindex!(cf[3],args[3],i...)
  f4 = getindex!(cf[4],args[4],i...)
  evaluate!(cgi,gi,f1,f2,f3,f4)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any,Any,Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  f3 = getindex!(cf[3],args[3],i...)
  f4 = getindex!(cf[4],args[4],i...)
  f5 = getindex!(cf[5],args[5],i...)
  evaluate!(cgi,gi,f1,f2,f3,f4,f5)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any,Any,Any,Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  f3 = getindex!(cf[3],args[3],i...)
  f4 = getindex!(cf[4],args[4],i...)
  f5 = getindex!(cf[5],args[5],i...)
  f6 = getindex!(cf[6],args[6],i...)
  evaluate!(cgi,gi,f1,f2,f3,f4,f5,f6)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any,Any,Any,Any,Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  f3 = getindex!(cf[3],args[3],i...)
  f4 = getindex!(cf[4],args[4],i...)
  f5 = getindex!(cf[5],args[5],i...)
  f6 = getindex!(cf[6],args[6],i...)
  f7 = getindex!(cf[7],args[7],i...)
  evaluate!(cgi,gi,f1,f2,f3,f4,f5,f6,f7)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any,Any,Any,Any,Any,Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  f3 = getindex!(cf[3],args[3],i...)
  f4 = getindex!(cf[4],args[4],i...)
  f5 = getindex!(cf[5],args[5],i...)
  f6 = getindex!(cf[6],args[6],i...)
  f7 = getindex!(cf[7],args[7],i...)
  f8 = getindex!(cf[8],args[8],i...)
  evaluate!(cgi,gi,f1,f2,f3,f4,f5,f6,f7,f8)
end

function _getindex_and_call!(cgi,gi,cf,args::Tuple{Any,Any,Any,Any,Any,Any,Any,Any,Any},i...)
  f1 = getindex!(cf[1],args[1],i...)
  f2 = getindex!(cf[2],args[2],i...)
  f3 = getindex!(cf[3],args[3],i...)
  f4 = getindex!(cf[4],args[4],i...)
  f5 = getindex!(cf[5],args[5],i...)
  f6 = getindex!(cf[6],args[6],i...)
  f7 = getindex!(cf[7],args[7],i...)
  f8 = getindex!(cf[8],args[8],i...)
  f9 = getindex!(cf[9],args[9],i...)
  evaluate!(cgi,gi,f1,f2,f3,f4,f5,f6,f7,f8,f9)
end
