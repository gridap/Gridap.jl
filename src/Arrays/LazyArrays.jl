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
@inline function lazy_map(k,f::AbstractArray...)
  s = common_size(f...)
  lazy_map(evaluate,Fill(k, s), f...)
end

@inline lazy_map(::typeof(evaluate),k::AbstractArray,f::AbstractArray...) = LazyArray(k,f...)

"""
    lazy_map(::Type{T},f,a::AbstractArray...) where T

Like [`lazy_map(f,a::AbstractArray...)`](@ref), but the user provides the element type
of the resulting array in order to circumvent type inference.
"""
@inline function lazy_map(k,T::Type,f::AbstractArray...)
  s = common_size(f...)
  lazy_map(evaluate,T,Fill(k, s), f...)
end

@inline lazy_map(::typeof(evaluate),T::Type,k::AbstractArray,f::AbstractArray...) = LazyArray(T,k,f...)

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
end

function LazyArray(g::AbstractArray{S}, f::AbstractArray...) where S
  isconcretetype(S) ? gi = testitem(g) : @notimplemented
  fi = map(testitem,f)
  T = return_type(gi, fi...)
  LazyArray(T, g, f...)
end

IndexStyle(::Type{<:LazyArray}) = IndexCartesian()

IndexStyle(::Type{<:LazyArray{G,T,1} where {G,T}}) = IndexLinear()

function array_cache(a::LazyArray)
  @notimplementedif ! all(map(isconcretetype, map(eltype, a.f)))
  if ! (eltype(a.g) <: Function)
    @notimplementedif ! isconcretetype(eltype(a.g))
  end
  gi = testitem(a.g)
  fi = map(testitem,a.f)
  cg = array_cache(a.g)
  cf = map(array_cache,a.f)
  cgi = return_cache(gi, fi...)
  cg, cgi, cf
end

@inline function getindex!(cache, a::LazyArray, i...)
  cg, cgi, cf = cache
  gi = getindex!(cg, a.g, i...)
  fi = map((ci,ai) -> getindex!(ci,ai,i...),cf,a.f)
  vi = evaluate!(cgi, gi, fi...)
  vi
end

function Base.getindex(a::LazyArray, i...)
  ca = array_cache(a)
  getindex!(ca, a, i...)
end

Base.size(a::LazyArray) = size(a.g)

# Particular implementations for Fill

function lazy_map(::typeof(evaluate),f::Fill, a::Fill...)
  ai = map(ai->ai.value,a)
  r = evaluate(f.value, ai...)
  s = common_size(f, a...)
  Fill(r, s)
end

function lazy_map(::typeof(evaluate),::Type{T}, f::Fill, a::Fill...) where T
  lazy_map(evaluate, f, a...)
end

# @santiagobadia : CompressedArray and Union{CompressedArray,Fill}
# To be done when starting Algebra part

#@fverdugo: I find the grad argument very confusing. It seems very specific for arrays of Maps/Fields
# whereas LazyArray is something more general.
# In fact, I don't believe we need this. It seems that it is not used in the tests, right?
# Perhaps, what you really need is something similar to test_array_of_fields of the old Gridap verison.
# Moreover, line
#  ax = lazy_map(a, x)
#  Seems outdated
function test_lazy_array(
  a::AbstractArray,
  x::AbstractArray,
  v::AbstractArray,
  cmp::Function=(==);
  grad=nothing)

  ax = lazy_map(a, x)
  test_array(ax, v, cmp)

  ca, cfi, cx = array_cache(a, x)
  t = true
  for i in 1:length(a)
    fi = getindex!(ca, a, i)
    xi = getindex!(cx, x, i)
    fxi = evaluate!(cfi, fi, xi)
    vi = v[i]
    ti = cmp(fxi, vi)
    t = t && ti
  end
  @test t

  if grad != nothing
    g = lazy_map(gradient, a)
    test_lazy_array(g, x, grad, cmp)
  end
end

function common_size(a::AbstractArray...)
  a1, = a
  @check all(map(ai->length(a1) == length(ai),a)) "Array sizes $(map(size,a)) are not compatible."
  if all( map(ai->size(a1) == size(ai),a) )
    size(a1)
  else
    (length(a1),)
  end
end

# struct ArrayWithCounter{T,N,A,C} <: AbstractArray{T,N}
#   array::A
#   counter::C
#   function ArrayWithCounter(a::AbstractArray{T,N}) where {T,N}
#     c = zeros(Int,size(a))
#     c[:] .= 0
#     new{T,N,typeof(a),typeof(c)}(a,c)
#   end
# end

# Base.size(a::ArrayWithCounter) = size(a.array)

# function Base.getindex(a::ArrayWithCounter,i::Integer...)
#   a.counter[i...] += 1
#   a.array[i...]
# end

# Base.IndexStyle(::Type{<:ArrayWithCounter{T,N,A}}) where {T,A,N} = IndexStyle(A)

# function reset_counter!(a::ArrayWithCounter)
#   a.counter[:] .= 0
# end

# @santiagobadia : Do we want the same names for both ?????
# @inline array_cache(a::LazyArray,i...) = return_cache(a,i...)
# @inline get_index!(c,a::LazyArray,i...) = evaluate!(c,a,i...)
