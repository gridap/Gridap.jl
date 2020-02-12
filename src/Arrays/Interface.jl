"""
    getindex!(cache,a::AbstractArray,i...)

Returns the item of the array `a` associated with index `i`
by (possibly) using the scratch data passed in the `cache` object.

It defaults to

    getindex!(cache,a::AbstractArray,i...) = a[i...]

# Examples

Iterating over an array using the `getindex!` function

```jldoctest
using Gridap.Arrays

a = collect(10:15)

cache = array_cache(a)
for i in eachindex(a)
  ai = getindex!(cache,a,i)
  println("\$i -> \$ai")
end

# output
1 -> 10
2 -> 11
3 -> 12
4 -> 13
5 -> 14
6 -> 15
```

"""
getindex!(cache,a::AbstractArray,i...) = a[i...]

"""
    array_cache(a::AbstractArray)

Returns a cache object to be used in the [`getindex!`](@ref) function.
It defaults to 

    array_cache(a::T) where T = nothing

for types `T` such that `uses_hash(T) == Val(false)`, and 

    function array_cache(a::T) where T
      hash = Dict{UInt,Any}()
      array_cache(hash,a)
    end

for types `T` such that `uses_hash(T) == Val(true)`, see the [`uses_hash`](@ref) function. In the later case, the
type `T` should implement the following signature:

    array_cache(hash::Dict,a::AbstractArray)

where we pass a dictionary (i.e., a hash table) in the first argument. This hash table can be used to test
if the object `a` has already build a cache and re-use it as follows

    id = objectid(a)
    if haskey(hash,id)
      cache = hash[id] # Reuse cache
    else
      cache = ... # Build a new cache depending on your needs
      hash[id] = cache # Register the cache in the hash table
    end

This mechanism is needed, e.g., to re-use intermediate results in complex lazy operation trees.
In multi-threading computations, a different hash table per thread has to be used in order
to avoid race conditions.
"""
function array_cache(a::AbstractArray)
  _default_array_cache(a,uses_hash(a))
end

function array_cache(hash,a::T) where T
  if uses_hash(T) == Val{true}()
    error("array_cache(::Dict,::$T) not defined")
  end
  array_cache(a)
end

function _default_array_cache(a,::Val{false})
  nothing
end

function _default_array_cache(a,::Val{true})
  hash = Dict{UInt,Any}()
  array_cache(hash,a)
end

"""
    uses_hash(::Type{<:AbstractArray})

This function is used to specify if the type `T` uses the
hash-based mechanism to reuse caches.  It should return
either `Val(true)` or `Val(false)`. It defaults to

    uses_hash(::Type{<:AbstractArray}) = Val(false)

Once this function is defined for the type `T` it can also
be called on instances of `T`.

"""
uses_hash(::Type{<:AbstractArray}) = Val(false)

uses_hash(::T) where T = uses_hash(T)

"""
$(TYPEDSIGNATURES)

Returns an arbitrary instance of `eltype(a)`. The default returned value is the first entry
in the array if `length(a)>0` and `testvalue(eltype(a))` if `length(a)==0`
See the [`testvalue`](@ref) function.

# Examples

```jldoctest
using Gridap.Arrays

a = collect(3:10)
ai = testitem(a)

b = Int[]
bi = testitem(b)

(ai, bi)

# output
(3, 0)

```

"""
function testitem(a::AbstractArray{T}) where T
  if length(a) >0
    first(a)
  else
    testvalue(T)
  end::T
end

function testitem(a::Fill)
  a.value
end

"""
    get_array(a::AbstractArray)

Returns `a`.
"""
function get_array(a::AbstractArray)
  a
end

# Test the interface

"""
    test_array(
      a::AbstractArray{T,N}, b::AbstractArray{S,N},cmp=(==)) where {T,S,N}

Checks if the entries in `a` and `b` are equal using the comparison function `cmp`.
It also stresses the new methods added to the `AbstractArray` interface interface.
"""
function test_array(
  a::AbstractArray{T,N}, b::AbstractArray{S,N},cmp=(==)) where {T,S,N}
  @test cmp(a,b)
  cache = array_cache(a)
  t = true
  for i in eachindex(a)
    bi = b[i]
    ai = getindex!(cache,a,i)
    t = t && cmp(bi,ai)
  end
  @test t
  t = true
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
    t = t && (typeof(ai) == eltype(a))
    t = t && (typeof(ai) == T)
  end
  @test t
  @test IndexStyle(a) == IndexStyle(b)
  true
end

# Some API

"""
    testitems(b::AbstractArray...) -> Tuple

Returns a tuple with the result of `testitem` applied to each of the
arrays in `b`.

# Examples

```jldoctest
using Gridap.Arrays

a = collect(3:10)
b = Int[]
c = Float64[]
d = ones(10)

testitems(a,b,c,d)

# output
(3, 0, 0.0, 1.0)

```
"""
function testitems(a::AbstractArray,b::AbstractArray...)
  va = testitem(a)
  vb = testitems(b...)
  (va,vb...)
end

function testitems(a::AbstractArray)
  va = testitem(a)
  (va,)
end

"""
    array_caches(a::AbstractArray...) -> Tuple

Returns a tuple with the cache of each array in `a`.
"""
function array_caches(a::AbstractArray,b::AbstractArray...)
  hash = Dict{UInt,Any}()
  array_caches(hash,a,b...)
end

function array_caches(hash::Dict,a::AbstractArray,b::AbstractArray...)
  ca = array_cache(hash,a)
  cb = array_caches(hash,b...)
  (ca,cb...)
end

function array_caches(hash::Dict,a::AbstractArray)
  ca = array_cache(hash,a)
  (ca,)
end

"""
    getitems!(c::Tuple,a::Tuple,i...) -> Tuple

Extracts the `i`-th entry of all arrays in the tuple `a` using the caches in the tuple
`c`. The results is a tuple containing each one of the extracted entries.

# Example

Iterating over three different arrays simultaneously using `getitems!`

```jldoctest
using Gridap.Arrays

a = collect(0:5)
b = collect(10:15)
c = collect(20:25)

caches = array_caches(a,b,c)
for i in eachindex(a)
   s = getitems!(caches,(a,b,c),i)
   println("\$i -> \$s")
end

# output
1 -> (0, 10, 20)
2 -> (1, 11, 21)
3 -> (2, 12, 22)
4 -> (3, 13, 23)
5 -> (4, 14, 24)
6 -> (5, 15, 25)
```

"""
@inline function getitems!(cf::Tuple,a::Tuple{Vararg{<:AbstractArray}},i...)
  _getitems!(cf,i,a...)
end

@inline function _getitems!(c,i,a,b...)
  ca,cb = _split(c...)
  ai = getindex!(ca,a,i...)
  bi = getitems!(cb,b,i...)
  (ai,bi...)
end

@inline function _getitems!(c,i,a)
  ca, = c
  ai = getindex!(ca,a,i...)
  (ai,)
end
