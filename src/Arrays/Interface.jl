"""
    getindex!(cache,a::AbstractArray,i...)

Returns the item of the array `a` associated with index `i`
by (possibly) using the scratch data passed in the `cache` object.

It defaults to

    getindex!(cache,a::AbstractArray,i...) = a[i...]

As for standard Julia arrays, the user needs to implement only one of the following signatures
depending on the `IndexStyle` of the array.

    getindex!(cache,a::AbstractArray,i::Integer)
    getindex!(cache,a::AbstractArray{T,N},i::Vararg{Integer,N}) where {T,N}

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
getindex!(cache,a::AbstractArray,i::CartesianIndex) = getindex!(cache,a,Tuple(i)...)
getindex!(cache,a::AbstractArray,i::Integer) = _getindex_1d!(IndexStyle(a),cache,a,i)
function getindex!(cache,a::AbstractArray{T,N},i::Vararg{Integer,N}) where {T,N}
  _getindex_nd!(IndexStyle(a),cache,a,CartesianIndex(i))
end
_getindex_1d!(s::IndexLinear,cache,a,i) = a[i]
_getindex_1d!(s::IndexCartesian,cache,a,i) = _getindex_nd!(s,cache,a,CartesianIndices(a)[i])
_getindex_nd!(s::IndexLinear,cache,a,i) = _getindex_1d!(s,cache,a,LinearIndices(a)[i])
_getindex_nd!(s::IndexCartesian,cache,a,i) = a[i]

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
if the object `a` has already built a cache and re-use it as follows

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
array_cache(a::AbstractArray) = _default_array_cache(a,uses_hash(a))
array_cache(hash::Dict,a::AbstractArray) = _default_array_cache(hash,a,uses_hash(a))
_default_array_cache(a,s::Val{true}) = array_cache(Dict{UInt,Any}(),a)
_default_array_cache(a,s::Val{false}) = nothing
_default_array_cache(hash::Dict,a,s::Val{false}) = array_cache(a)
_default_array_cache(hash::Dict,a,s::Val{true}) = @abstractmethod

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
  #@check isconcretetype(T) "This array is type-instable"
  if length(a) >0
    first(a)
  else
    testvalue(T)
  end::T
end

function testitem(a::Fill)
  a.value
end

function testitem(a::Number)
  a
end

"""
    testvalue(::Type{T}) where T

Returns an arbitrary instance of type `T`. It defaults to `zero(T)` for
non-array types and to an empty array for array types.
It can be overloaded for new types `T` if `zero(T)` does not makes sense.
This function is used to compute  [`testitem`](@ref) for 0-length arrays.
"""
function testvalue end

testvalue(::Type{T}) where T = zero(T)
testvalue(v) = testvalue(typeof(v))

function testvalue(::Type{T}) where T<:AbstractArray{E,N} where {E,N}
   similar(T,tfill(0,Val(N))...)
end

# When the jacobian of a residual is obtained through automatic differentiation,
# the return type is BlockArray{<:SubArray} and the behaviour of testvalue
# does not allow broadcasting operations between BlockArray{<:AbstractMatrix}
# and BlockArray{<:SubArray}. This function returns a matrix of size a
# P-dimensional array where each dimension has length 0, i.e., (0, ..., 0).
function testvalue(::Type{<:SubArray{T,P,AT}}) where {T,P,AT}
  a = testvalue(AT)
  return SubArray(a, ntuple(_ -> 0:-1, P))
end

function testvalue(::Type{T}) where T<:Transpose{E,A} where {E,A}
  a = testvalue(A)
  Transpose(a)
end

testvalue(::Type{Base.OneTo{T}}) where T = Base.OneTo(zero(T))

testvalue(::Type{Base.UnitRange{T}}) where T = UnitRange(one(T),zero(T))

function testvalue(::Type{T}) where T<:Fill{E,N,A} where {E,N,A}
  Fill(zero(E),testvalue(A))
end

function testvalue(::Type{<:Tuple})
  @notimplemented "testvalue on Tuple type only implemented up to 8 tuple elements"
end

#@fverdugo: use meta-programming here
function testvalue(::Type{Tuple{T1,T2,T3,T4,T5,T6,T7,T8}}) where {T1,T2,T3,T4,T5,T6,T7,T8}
  (testvalue(T1),testvalue(T2),testvalue(T3),testvalue(T4),testvalue(T5),testvalue(T6),testvalue(T7),testvalue(T8))
end

function testvalue(::Type{Tuple{T1,T2,T3,T4,T5,T6,T7}}) where {T1,T2,T3,T4,T5,T6,T7}
  (testvalue(T1),testvalue(T2),testvalue(T3),testvalue(T4),testvalue(T5),testvalue(T6),testvalue(T7))
end

function testvalue(::Type{Tuple{T1,T2,T3,T4,T5,T6}}) where {T1,T2,T3,T4,T5,T6}
  (testvalue(T1),testvalue(T2),testvalue(T3),testvalue(T4),testvalue(T5),testvalue(T6))
end

function testvalue(::Type{Tuple{T1,T2,T3,T4,T5}}) where {T1,T2,T3,T4,T5}
  (testvalue(T1),testvalue(T2),testvalue(T3),testvalue(T4),testvalue(T5))
end

function testvalue(::Type{Tuple{T1,T2,T3,T4}}) where {T1,T2,T3,T4}
  (testvalue(T1),testvalue(T2),testvalue(T3),testvalue(T4))
end

function testvalue(::Type{Tuple{T1,T2,T3}}) where {T1,T2,T3}
  (testvalue(T1),testvalue(T2),testvalue(T3))
end

function testvalue(::Type{Tuple{T1,T2}}) where {T1,T2}
  (testvalue(T1),testvalue(T2))
end

function testvalue(::Type{Tuple{T1}}) where {T1}
  (testvalue(T1),)
end

function testvalue(::Type{Tuple{}})
  ()
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
It also stresses the new methods added to the `AbstractArray` interface.
"""
function test_array(
  a::AbstractArray{T,N}, b::AbstractArray{S,N},cmp=(==)) where {T,S,N}

  function _test_loop(indices)
    cache = array_cache(a)
    t = true
    for i in indices
      bi = b[i]
      ai = getindex!(cache,a,i)
      t = t && cmp(bi,ai)
    end
    @test t
  end

  @test cmp(a,b)
  _test_loop(eachindex(a))
  _test_loop(LinearIndices(a))
  _test_loop(CartesianIndices(a))
  cache = array_cache(a)
  t = true
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
    t = t && (typeof(ai) <: eltype(a))
    t = t && (typeof(ai) <: T)
  end
  @test t
  @test isa(testitem(a),eltype(a))
  if length(a) > 0
    @test testitem(a) == first(a)
  end
  true
end
