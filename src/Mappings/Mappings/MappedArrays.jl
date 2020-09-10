"""
    apply(f,a::AbstractArray...) -> AbstractArray

Applies the mapping `f` to the entries of the arrays in `a` (see the definition of [`mapping`](@ref)).

The resulting array `r` is such that `r[i]` equals to `apply_mapping(f,ai...)` where `ai`
is the tuple containing the `i`-th entry of the arrays in `a` (see function
[`apply_mapping`](@ref) for more details).
In other words, the resulting array is numerically equivalent to:

    map( (x...)->apply_mapping(f,x...), a...)


See the [`apply_mapping`](@ref) function for details.

# Examples

Using a function as mapping

```jldoctest
using Gridap.Arrays

a = collect(0:5)
b = collect(10:15)

c = apply(+,a,b)

println(c)

# output
[10, 12, 14, 16, 18, 20]
```

Using a user-defined mapping

```jldoctest
using Gridap.Arrays
import Gridap.Arrays: apply_mapping!

a = collect(0:5)
b = collect(10:15)

struct MySum <: mapping end

apply_mapping!(cache,::MySum,x,y) = x + y

k = MySum()

c = apply(k,a,b)

println(c)

# output
[10, 12, 14, 16, 18, 20]
```

"""

apply_mapping(f::Function,a...) = _apply_broadcast(f,a...)

apply_mapping(f::Mapping,a...) = _apply_broadcast(f,a...)

function _apply_broadcast(f,a::AbstractArray...)
    s = common_size(a...)
    apply_mapping(Fill(f, s...), a...)
end

"""
    apply(::Type{T},f,a::AbstractArray...) where T

Like [`apply(f,a::AbstractArray...)`](@ref), but the user provides the element type
of the resulting array in order to circumvent type inference.
"""
function apply_mapping(::Type{T}, f, a::AbstractArray...) where T
    s = common_size(a...)
    apply_mapping(T, Fill(f, s...), a...)
end

"""
    apply(f::AbstractArray,a::AbstractArray...) -> AbstractArray
Applies the mappings in the array of mappings `f` to the entries in the arrays in `a`.

The resulting array has the same entries as the one obtained with:

    map( apply_mapping, f, a...)

See the [`apply_mapping`](@ref) function for details.

# Example

"Evaluating" an array of functions

```jldoctest
using Gridap.Arrays

f = [+,-,max,min]
a = [1,2,3,4]
b = [4,3,2,1]

c = apply(f,a,b)

println(c)

# output
[5, -1, 3, 1]
```
"""
function apply_mapping(f::AbstractArray, a::AbstractArray...)
    MappedArray(f, a...)
end

"""
    apply(::Type{T},f::AbstractArray,a::AbstractArray...) where T

Like [`apply(f::AbstractArray,a::AbstractArray...)`](@ref), but the user provides the element type
of the resulting array in order to circumvent type inference.
"""
function apply_mapping(::Type{T}, f::AbstractArray, a::AbstractArray...) where T
    MappedArray(T, f, a...)
end

# @santiagobadia : Do we really need this?
function apply_mappings(fs::Tuple, x::AbstractArray)
    return Tuple([apply_mapping(f, ax) for f in fs])
end

# Helpers

struct MappedArray{T,N,F,G} <: AbstractArray{T,N}
    g::G
    f::F
    function MappedArray(::Type{T}, g::AbstractArray, f::AbstractArray...) where T
        G = typeof(g)
        F = typeof(f)
        f1, = f
        new{T,ndims(f1),F,G}(g, f)
    end
end

function MappedArray(g::AbstractArray, f::AbstractArray...)
    gi = testitem(g) # Assumes that all mappings return the same type
    fi = testitems(f...)
    T = typeof(testitem(gi, fi...))
    MappedArray(T, g, f...)
end

function uses_hash(::Type{<:MappedArray})
    Val(true)
end

function array_cache(hash::Dict, a::MappedArray)
    id = objectid(a)
    _cache = _array_cache(hash, a)
    if haskey(hash, id)
        cache = _get_stored_cache(_cache, hash, id)
    else
        cache = _cache
        hash[id] = cache
    end
    cache
end

# TODO
@noinline function _get_stored_cache(cache::T, hash, id) where T
    hash[id]
end

function _array_cache(hash, a::MappedArray)
    cg = array_cache(hash, a.g)
    gi = testitem(a.g)
    fi = testitems(a.f...)
    @notimplementedif ! all(map(isconcretetype, map(eltype, a.f)))
    if ! (eltype(a.g) <: Function)
        @notimplementedif ! isconcretetype(eltype(a.g))
    end
    cf = array_caches(hash, a.f...)
    cgi = return_cache(gi, fi...)
    ai = testitem!(cgi, gi, fi...)
    i = -testitem(eachindex(a))
    e = Evaluation((i,), ai)
    c = (cg, cgi, cf)
    (c, e)
end

function testitem(a::MappedArray)
    cg = array_cache(a.g)
    gi = testitem(a.g)
    fi = testitems(a.f...)
    testitem(gi, fi...)
end

function getindex!(cache, a::MappedArray, i::Integer...)
    li = LinearIndices(a)
    getindex!(cache, a, li[i...])
end

function getindex!(cache, a::MappedArray, i::Integer)
    _cached_getindex!(cache, a, (i,))
end

function getindex!(cache, a::MappedArray, i::CartesianIndex)
    _cached_getindex!(cache, a, Tuple(i))
end

function _cached_getindex!(cache, a::MappedArray, i::Tuple)
    c, e = cache
    v = e.fx
    if e.x != i
        v = _getindex!(c, a, i...)
        e.x = i
        e.fx = v
    end
    v
end

function _getindex!(cache, a::MappedArray, i...)
    cg, cgi, cf = cache
    gi = getindex!(cg, a.g, i...)
    fi = getitems!(cf, a.f, i...)
    vi = evaluate!(cgi, gi, fi...)
    vi
end

function Base.getindex(a::MappedArray, i...)
    ca = array_cache(a)
    getindex!(ca, a, i...)
end

function Base.IndexStyle(
  ::Type{MappedArray{T,N,F,G}}) where {T,N,F,G}
    common_index_style(F)
end

function Base.size(a::MappedArray)
    f, = a.f
    size(f)
end

function common_size(a::AbstractArray...)
    a1, = a
    c = all([length(a1) == length(ai) for ai in a])
    if !c
        error("Array sizes $(map(size, a)) are not compatible.")
    end
    l = length(a1)
    (l,)
end

# TODO not sure what to do with shape and index-style
function common_index_style(::Type{<:Tuple})
    IndexLinear()
end

mutable struct Evaluation{X,F} <: GridapType
    x::X
    fx::F
    function Evaluation(x::X, fx::F) where {X,F}
        new{X,F}(x, fx)
    end
end

# Particular implementations for Fill

function apply_mapping(f::Fill, a::Fill...)
    ai = _getvalues(a...)
    r = evaluate(f.value, ai...)
    s = common_size(f, a...)
    Fill(r, s)
end

function apply_mapping(::Type{T}, f::Fill, a::Fill...) where T
    apply_mapping(f, a...)
end

function _getvalues(a::Fill, b::Fill...)
    ai = a.value
    bi = _getvalues(b...)
    (ai, bi...)
end

function _getvalues(a::Fill)
    ai = a.value
    (ai,)
end

function Base.sum(a::MappedArray)
    cache = array_cache(a)
    _sum_array_cache(cache, a)
end

function _sum_array_cache(cache, a)
    r = zero(eltype(a))
    for i in eachindex(a)
        ai = getindex!(cache, a, i)
        r += ai
    end
    r
end

# TODO Think about iteration and sub-iteration

# Helper type to test caching of intermediate results

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

"""
    field_array_cache(a::AbstractArray,x::AbstractArray) -> Tuple

    Returns the caches needed to perform the following iteration

        ca, cfi, cx = field_array_cache(a,x)

        for i in length(a)
          fi = getindex!(ca,a,i)
          xi = getindex!(cx,x,i)
          fxi = evaluate!(cfi,fi,xi)
        end
"""
function return_mapping_array_cache(a::AbstractArray, x::AbstractArray)
  ca = array_cache(a)
  fi = testitem(a)
  xi = testitem(x)
  cfi = return_cache(fi, xi)
  cx = array_cache(x)
  (ca, cfi, cx)
end

function test_mapped_array(
  a::AbstractArray,
  x::AbstractArray,
  v::AbstractArray,
  cmp::Function=(==);
  grad=nothing)

  ax = apply_mapping(a, x)
  test_array(ax, v, cmp)

  ca, cfi, cx = return_mapping_array_cache(a, x)
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
    g = apply_mapping(gradient, a)
    test_mapped_array(g, x, grad, cmp)
  end

end
