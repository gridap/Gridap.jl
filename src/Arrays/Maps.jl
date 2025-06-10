"""
Abstract type representing a function (mapping) that provides a cache and an in-place
evaluation for performance. This is the type to be used in the [`lazy_map`](@ref) function.

Derived types must implement the following method:

- [`evaluate!(cache,k,x...)`](@ref)

and optionally these ones:

- [`return_cache(k,x...)`](@ref)
- [`return_type(k,x...)`](@ref)

The mapping interface can be tested with the [`test_map`](@ref) function.

Note that most of the functionality implemented in terms of this interface
relies in duck typing. That is, it is not strictly needed to work with types
that inherit from `Map`. This is specially useful in order to accommodate
existing types into this framework without the need to implement a wrapper type
that inherits from `Map`. For instance, a default implementation is available
for `Function` objects.  However, we recommend that new types inherit from `Map`.

"""
abstract type Map <: GridapType end

"""
    return_cache(f,x...)

Returns the `cache` needed to lazy_map mapping `f` with arguments
of the same type as the objects in `x`.
This function returns `nothing` by default, i.e., no cache.
"""
return_cache(f,x...) = nothing

"""
    evaluate!(cache,f,x...)

Applies the mapping `f` at the arguments `x...` using
the scratch data provided in the given `cache` object. The `cache` object
is built with the [`return_cache`](@ref) function using arguments of the same type as in `x`.
In general, the returned value `y` can share some part of its state with the `cache` object.
If the result of two or more calls to this function need to be accessed simultaneously
(e.g., in multi-threading), create and use several `cache` objects (e.g., one cache
per thread).
"""
@noinline function evaluate!(cache,f,x...)
  @abstractmethod """
  Method

      Gridap.Arrays.evaluate!(cache,f,x...)

  has not been defined for the requested types (see stack trace).
  """
end

"""
    return_type(f,x...)

Returns the type of the result of calling mapping `f` with
arguments of the types of the objects `x`.
"""
return_type(f,x...) = typeof(return_value(f,x...))

"""
    return_value(f,x...)

Return a variable of the type of the image fx=`f`(`x`...) (possibly fx itself).
"""
return_value(f,x...) = evaluate(f,testargs(f,x...)...)

"""
    testargs(f,x...)

The default implementation of this function is `testargs(f,x...) = x`.
One can overload it in order to use `lazy_map` with 0-length array
and maps with non-trivial domains.
"""
testargs(f,x...) = x

"""
    evaluate(f,x...)

evaluates the mapping `f` at the arguments in `x` by creating a temporary cache
internally. This functions is equivalent to
```jl
cache = return_cache(f,x...)
evaluate!(cache,f,x...)
```
"""
function evaluate(f,x...)
  c = return_cache(f,x...)
  y = evaluate!(c,f,x...)
end

(m::Map)(x...) = evaluate(m,x...)

# Default implementation for Function
evaluate!(cache,f::Function,x...) = f(x...)

# Default implementation for constructors
evaluate!(cache,::Type{f},x...) where f = f(x...)


# Testing the interface
"""
    test_map(y,f,x...;cmp=(==))

Function used to test if the mapping `f` has been
implemented correctly. `f` is a `Map` sub-type, `x` is a tuple in the domain of the
mapping and `y` is the expected result. Function `cmp` is used to compare
the computed result with the expected one. The checks are done with the `@test`
macro.
"""
function test_map(y,f,x...;cmp=(==))
  z = evaluate(f,x...)
  @test cmp(z,y)
  @test typeof(z) == return_type(f,x...)
  cache = return_cache(f,x...)
  z = evaluate!(cache,f,x...)
  @test cmp(z,y)
  z = evaluate!(cache,f,x...)
  @test cmp(z,y)
  true
end

# Broadcast Functions

"""
    Broadcasting(f)

Returns a mapping that represents the "broadcasted" version of the
function `f`.

# Example

```jldoctest
using Gridap.Arrays

a = [3,2]
b = [2,1]

bm = Broadcasting(+)

c = evaluate(bm,a,b)

println(c)

# output
[5, 3]
```
"""
struct Broadcasting{F} <: Map
  f::F
end

return_cache(f::Broadcasting,x...) = nothing

evaluate!(cache,f::Broadcasting,x...) = broadcast(f.f,x...)

function return_value(f::Broadcasting,x...)
  broadcast( (y...) -> f.f(testargs(f.f,y...)...), x... )
end

function evaluate!(cache,f::Broadcasting,x::Union{Number,AbstractArray{<:Number}}...)
  r = _prepare_cache!(cache,x...)
  a = r.array
  broadcast!(f.f,a,x...)
  a
end

function evaluate!(cache,f::Broadcasting,x::AbstractArray{<:Number})
  setsize!(cache,size(x))
  a = cache.array
  @inbounds for i in eachindex(x)
    a[i] = f.f(x[i])
  end
  a
end

function evaluate!(cache,f::Broadcasting,args::Number...)
  f.f(args...)
end

function return_value(f::Broadcasting,x::Number...)
  return_value(f.f,x...)
end

function return_cache(f::Broadcasting,x::Number...)
  nothing
end

function return_value(f::Broadcasting,x::Union{Number,AbstractArray{<:Number}}...)
  s = map(_size_zero,x)
  bs = Base.Broadcast.broadcast_shape(s...)
  T = return_type(f.f,map(testitem,x)...)
  r = fill(testvalue(T),bs)
  r
end

function _size_zero(a)
  s = size(a)
  if length(a) == 0
    r = map(i-> (i==0 ? 1 : i) ,s)
  else
    r = s
  end
  r
end
_size_zero(a::Number) = (1,)

function return_cache(f::Broadcasting,x::Union{Number,AbstractArray{<:Number}}...)
  s = map(_size,x)
  bs = Base.Broadcast.broadcast_shape(s...)
  T = return_type(f.f,map(testitem,x)...)
  r = fill(testvalue(T),bs)
  cache = CachedArray(r)
  _prepare_cache!(cache,x...)
  cache
end

function _prepare_cache!(c,x...)
  s = map(_size,x)
  bs = Base.Broadcast.broadcast_shape(s...)
  if bs != size(c)
    setsize!(c,bs)
  end
  c
end

_size(a) = size(a)
_size(a::Number) = (1,)

"""
    OperationMap(f,args)

Returns a mapping that represents the result of applying the function `f`
to the arguments in the tuple `args`.
That is, `OperationMap(f,args)(x...)` is formally defined as
`f(map(a->a(x...),args)...)`
"""
struct OperationMap{K,L} <: Map
  k::K
  l::L
  function OperationMap(k,l)
    new{typeof(k),typeof(l)}(k,l)
  end
end

function return_cache(c::OperationMap,x...)
  cl = map(fi -> return_cache(fi,x...),c.l)
  lx = map(fi -> return_value(fi,x...),c.l)
  ck = return_cache(c.k,lx...)
  ck, cl
end

function evaluate!(cache,c::OperationMap,x...)
  ck, cf = cache
  lx = map((ci,fi) -> evaluate!(ci,fi,x...),cf,c.l)
  evaluate!(ck,c.k,lx...)
end

# Operations

"""
    Operation(op)

Returns the map that results after applying an operation `f` over a set of map(s) `args`.
That is `Operation(f)(args)(x...)` is formally defined as
`f(map(a->a(x...),args)...)`.

# Example

```jldoctest
using Gridap.Arrays

fa(x) = x.*x
fb(x) = sqrt.(x)

x = collect(0:5)

fab = Operation(fa)(fb)
c = evaluate(fab,x)

println(c)

# output
[0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
```
"""
struct Operation{T} <: Map
  op::T
end

evaluate!(cache,op::Operation,args::Map...) = OperationMap(op.op,args)
evaluate!(cache,op::Operation,args::Function...) = OperationMap(op.op,args)

"""
"""
function inverse_map(f)
  @unreachable """\n
  Function inverse_map is not implemented yet for objects of type $(typeof(f))
  """
end

"""
    struct InverseMap{F} <: Map

Map for the inverse of the `Function` or [`Map`](@ref) `F`.
"""
struct InverseMap{F} <: Map
  original::F
end

function evaluate!(cache,k::InverseMap,args...)
  @notimplemented """\n
  The inverse evaluation is not implemented yet for maps of type $(typeof(k.original))
  """
end

inverse_map(k::Map) = InverseMap(k)
inverse_map(k::InverseMap) = k.original
