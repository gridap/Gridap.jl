"""
Abstract type representing the operations to be used in the [`apply`](@ref) function.

Derived types must implement the following method:

- [`apply_Mapping!(cache,k,x...)`](@ref)

and optionally these ones:

- [`Mapping_cache(k,x...)`](@ref)
- [`Mapping_return_type(k,x...)`](@ref)

The Mapping interface can be tested with the [`test_Mapping`](@ref) function.

Note that most of the functionality implemented in terms of this interface
relies in duck typing. That is, it is not strictly needed to work with types
that inherit from `Mapping`. This is specially useful in order to accommodate
existing types into this framework without the need to implement a wrapper type
that inherits from `Mapping`. For instance, a default implementation is available
for `Function` objects.  However, we recommend that new types inherit from `Mapping`.

"""
abstract type Mapping <: GridapType end

# const MappingOrFunction = Union{Mapping,Function}

"""
    mapping_return_type(f,x...)

Returns the type of the result of calling mapping `f` with
arguments of the types of the objects `x`.

It defaults to `typeof(Mapping_testitem(f,x...))`
"""
return_type(f,x...) = typeof(testitem(f,x...))

"""
    mapping_cache(f,x...)

Returns the `cache` needed to apply mapping `f` with arguments
of the same type as the objects in `x`.
This function returns `nothing` by default.
"""
return_cache(f,x...) = nothing

"""
    apply_mapping!(cache,f,x...)

Applies the mapping `f` at the arguments `x...` using
the scratch data provided in the given `cache` object. The `cache` object
is built with the [`mapping_cache`](@ref) function using arguments of the same type as in `x`.
In general, the returned value `y` can share some part of its state with the `cache` object.
If the result of two or more invocations of this function need to be accessed simultaneously
(e.g., in multi-threading), create and use various `cache` objects (e.g., one cache
per thread).
"""
evaluate!(cache,f,x...) = @abstractmethod

"""
    apply_mapping(f,x...)

apply the mapping `f` at the arguments in `x` by creating a temporary cache
internally. This functions is equivalent to
```jl
cache = mapping_cache(f,x...)
apply_mapping!(cache,f,x...)
```
"""
function evaluate(f,x...)
  c = return_cache(f,x...)
  y = evaluate!(c,f,x...)
  y
end

# @santiagobadia: Duck typing not ok here I guess
function testitem(k::Mapping,x...)
  cache = return_cache(k,x...)
  testitem!(cache,k,x...)
end

@inline function testitem!(cache,k::Mapping,x...)
  evaluate!(cache,k,x...)
end

# Differentiation

# @santiagobadia : We could consider autodiff here
derivative(k::Mapping) = @notimplemented

gradient(k::Mapping) = transpose(derivative(k::Mapping))

const ∇ = gradient

# Testing the interface

"""
    test_mapping(f,x::Tuple,y,cmp=(==))

Function used to test if the mapping `f` has been
implemented correctly. `f` is a mapping object, `x` is a tuple containing the arguments
of the mapping, and `y` is the expected result. Function `cmp` is used to compare
the computed result with the expected one. The checks are done with the `@test`
macro.
"""
function test_mapping(f,x::Tuple,y,cmp=(==))
  z = evaluate(f,x...)
  @test cmp(z,y)
  @test typeof(z) == return_type(f,x...)
  cache = return_cache(f,x...)
  z = evaluate!(cache,f,x...)
  @test cmp(z,y)
  z = evaluate!(cache,f,x...)
  @test cmp(z,y)
  z = testitem!(cache,f,x...)
  @test cmp(typeof(z),typeof(y))
end

# Work with several mappings at once

"""
    mapping_caches(fs::Tuple,x...) -> Tuple

Returns a tuple with the cache corresponding to each mapping in `fs`
for the arguments `x...`.
"""
function return_caches(fs::Tuple,x...)
  _mapping_caches(x,fs...)
end

function _mapping_caches(x::Tuple,a,b...)
  ca = return_cache(a,x...)
  cb = return_caches(b,x...)
  (ca,cb...)
end

function _mapping_caches(x::Tuple,a)
  ca = return_cache(a,x...)
  (ca,)
end

"""
    apply_mappings!(caches::Tuple,fs::Tuple,x...) -> Tuple

Applies the mappings in the tuple `fs` at the arguments `x...`
by using the corresponding cache objects in the tuple `caches`.
The result is also a tuple containing the result for each mapping in `fs`.
"""
@inline function evaluate!(cfs::Tuple,f::Tuple,x...)
  _evaluate_mappings!(cfs,x,f...)
end

@inline function _evaluate_mappings!(cfs,x,f1,f...)
  cf1, cf = _split(cfs...)
  f1x = evaluate!(cf1,f1,x...)
  fx = evaluate!(cf,f,x...)
  (f1x,fx...)
end

@inline function _evaluate_mappings!(cfs,x,f1)
  cf1, = cfs
  f1x = evaluate!(cf1,f1,x...)
  (f1x,)
end

@inline function _split(a,b...)
  (a,b)
end

"""
    mapping_return_types(f::Tuple,x...) -> Tuple

Computes the return types of the mappings in `f` when called
with arguments `x`.
"""
function return_types(f::Tuple,x...)
  _mapping_return_types(x,f...)
end

function _mapping_return_types(x::Tuple,a,b...)
  Ta = return_type(a,x...)
  Tb = return_types(b,x...)
  (Ta,Tb...)
end

function _mapping_return_types(x::Tuple,a)
  Ta = return_type(a,x...)
  (Ta,)
end

# Function implementation

evaluate!(cache,f::Function,x...) = f(x...)

return_cache(f::Function,x...) = nothing

evaluate(f::Function,x...) = f(x...)

evaluate(T::Type,f::Function,x...) = f(x...)

testitem(f::Function,x...) = f(x...)

testitem!(cache,f::Function,x...) = f(x...)

# Number or Array implementation

const NumberOrArray = Union{Number,AbstractArray{<:Number}}

evaluate!(cache,f::NumberOrArray,x...) = f

return_cache(f::NumberOrArray,x...) = nothing

return_type(f::NumberOrArray,x...) = typeof(f)

evaluate(f::NumberOrArray,x...) = f

evaluate(T::Type,f::NumberOrArray,x...) = f
