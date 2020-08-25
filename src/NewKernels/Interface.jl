
"""
Abstract type representing the operations to be used in the [`apply`](@ref) function.

Derived types must implement the following method:

- [`apply_kernel!(cache,k,x...)`](@ref)

and optionally these ones:

- [`kernel_cache(k,x...)`](@ref)
- [`kernel_return_type(k,x...)`](@ref)

The kernel interface can be tested with the [`test_kernel`](@ref) function.

Note that most of the functionality implemented in terms of this interface
relies in duck typing. That is, it is not strictly needed to work with types
that inherit from `Kernel`. This is specially useful in order to accommodate
existing types into this framework without the need to implement a wrapper type
that inherits from `Kernel`. For instance, a default implementation is available
for `Function` objects.  However, we recommend that new types inherit from `Kernel`.

"""
abstract type NewKernel <: GridapType end

"""
    kernel_return_type(f,x...)

Returns the type of the result of calling kernel `f` with
arguments of the types of the objects `x`.

It defaults to `typeof(kernel_testitem(f,x...))`
"""
return_type(f::NewKernel,x...) = typeof(testitem(f,x...))

"""
    kernel_cache(f,x...)

Returns the `cache` needed to apply kernel `f` with arguments
of the same type as the objects in `x`.
This function returns `nothing` by default.
"""
return_cache(f::NewKernel,x...) = nothing

"""
    apply_kernel!(cache,f,x...)

Applies the kernel `f` at the arguments `x...` using
the scratch data provided in the given `cache` object. The `cache` object
is built with the [`kernel_cache`](@ref) function using arguments of the same type as in `x`.
In general, the returned value `y` can share some part of its state with the `cache` object.
If the result of two or more invocations of this function need to be accessed simultaneously
(e.g., in multi-threading), create and use various `cache` objects (e.g., one cache
per thread).
"""
evaluate!(cache,f::NewKernel,x...) = @abstractmethod


# Testing the interface

"""
    test_kernel(f,x::Tuple,y,cmp=(==))

Function used to test if the kernel `f` has been
implemented correctly. `f` is a kernel object, `x` is a tuple containing the arguments
of the kernel, and `y` is the expected result. Function `cmp` is used to compare
the computed result with the expected one. The checks are done with the `@test`
macro.
"""
function test_kernel(f,x::Tuple,y,cmp=(==))
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


# Functions on kernel objects

"""
    apply_kernel(f,x...)

apply the kernel `f` at the arguments in `x` by creating a temporary cache
internally. This functions is equivalent to
```jl
cache = kernel_cache(f,x...)
apply_kernel!(cache,f,x...)
```
"""
function evaluate(f::NewKernel,x...)
  c = return_cache(f,x...)
  y = evaluate!(c,f,x...)
  y
end

# Work with several kernels at once

"""
    kernel_caches(fs::Tuple,x...) -> Tuple

Returns a tuple with the cache corresponding to each kernel in `fs`
for the arguments `x...`.
"""
function return_caches(fs::NTuple{N,<:NewKernel},x...) where N
  _kernel_caches(x,fs...)
end

function _kernel_caches(x::Tuple,a,b...)
  ca = return_cache(a,x...)
  cb = return_caches(b,x...)
  (ca,cb...)
end

function _kernel_caches(x::Tuple,a)
  ca = return_cache(a,x...)
  (ca,)
end

"""
    apply_kernels!(caches::Tuple,fs::Tuple,x...) -> Tuple

Applies the kernels in the tuple `fs` at the arguments `x...`
by using the corresponding cache objects in the tuple `caches`.
The result is also a tuple containing the result for each kernel in `fs`.
"""
@inline function evaluate!(cfs::Tuple,f::NTuple{N,<:NewKernel},x...) where N
  _evaluate_kernels!(cfs,x,f...)
end

@inline function _evaluate_kernels!(cfs,x,f1,f...)
  cf1, cf = _split(cfs...)
  f1x = evaluate!(cf1,f1,x...)
  fx = evaluate!(cf,f,x...)
  (f1x,fx...)
end

@inline function _evaluate_kernels!(cfs,x,f1)
  cf1, = cfs
  f1x = evaluate!(cf1,f1,x...)
  (f1x,)
end

@inline function _split(a,b...)
  (a,b)
end

"""
    kernel_return_types(f::Tuple,x...) -> Tuple

Computes the return types of the kernels in `f` when called
with arguments `x`.
"""
function return_types(f::Tuple,x...)
  _kernel_return_types(x,f...)
end

function _kernel_return_types(x::Tuple,a,b...)
  Ta = return_type(a,x...)
  Tb = return_types(b,x...)
  (Ta,Tb...)
end

function _kernel_return_types(x::Tuple,a)
  Ta = return_type(a,x...)
  (Ta,)
end

function testitem(k::NewKernel,x...)
  cache = return_cache(k,x...)
  testitem!(cache,k,x...)
end

@inline function testitem!(cache,k::NewKernel,x...)
  evaluate!(cache,k,x...)
end

# Some particular cases

const NumberOrArray = Union{Number,AbstractArray}
