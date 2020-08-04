# Define kernel interface

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
abstract type Kernel <: GridapType end

"""
    kernel_return_type(f,x...)

Returns the type of the result of calling kernel `f` with
arguments of the types of the objects `x`.

It defaults to `typeof(kernel_testitem(f,x...))`
"""
function kernel_return_type(f,x...)
  typeof(kernel_testitem(f,x...))
end

"""
    kernel_cache(f,x...)

Returns the `cache` needed to apply kernel `f` with arguments
of the same type as the objects in `x`.
This function returns `nothing` by default.
"""
kernel_cache(f,x...) = nothing

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
function apply_kernel!(cache,f,x...)
  @abstractmethod
end

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
  z = apply_kernel(f,x...)
  @test cmp(z,y)
  @test typeof(z) == kernel_return_type(f,x...)
  cache = kernel_cache(f,x...)
  z = apply_kernel!(cache,f,x...)
  @test cmp(z,y)
  z = apply_kernel!(cache,f,x...)
  @test cmp(z,y)
  z = kernel_testitem!(cache,f,x...)
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
function apply_kernel(f,x...)
  cache = kernel_cache(f,x...)
  y = apply_kernel!(cache,f,x...)
  y
end

# Work with several kernels at once

"""
    kernel_caches(fs::Tuple,x...) -> Tuple

Returns a tuple with the cache corresponding to each kernel in `fs`
for the arguments `x...`.
"""
function kernel_caches(fs::Tuple,x...)
  _kernel_caches(x,fs...)
end

function _kernel_caches(x::Tuple,a,b...)
  ca = kernel_cache(a,x...)
  cb = kernel_caches(b,x...)
  (ca,cb...)
end

function _kernel_caches(x::Tuple,a)
  ca = kernel_cache(a,x...)
  (ca,)
end

"""
    apply_kernels!(caches::Tuple,fs::Tuple,x...) -> Tuple

Applies the kernels in the tuple `fs` at the arguments `x...`
by using the corresponding cache objects in the tuple `caches`.
The result is also a tuple containing the result for each kernel in `fs`.
"""
@inline function apply_kernels!(cfs::Tuple,f::Tuple,x...)
  _apply_kernels!(cfs,x,f...)
end

@inline function _apply_kernels!(cfs,x,f1,f...)
  cf1, cf = _split(cfs...)
  f1x = apply_kernel!(cf1,f1,x...)
  fx = apply_kernels!(cf,f,x...)
  (f1x,fx...)
end

@inline function _apply_kernels!(cfs,x,f1)
  cf1, = cfs
  f1x = apply_kernel!(cf1,f1,x...)
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
function kernel_return_types(f::Tuple,x...)
  _kernel_return_types(x,f...)
end

function _kernel_return_types(x::Tuple,a,b...)
  Ta = kernel_return_type(a,x...)
  Tb = kernel_return_types(b,x...)
  (Ta,Tb...)
end

function _kernel_return_types(x::Tuple,a)
  Ta = kernel_return_type(a,x...)
  (Ta,)
end

function kernel_testitem(k,x...)
  cache = kernel_cache(k,x...)
  kernel_testitem!(cache,k,x...)
end

@inline function kernel_testitem!(cache,k,x...)
  apply_kernel!(cache,k,x...)
end

# Include some well-known types in this interface

function kernel_return_type(f::Function,x...)
  Ts = map(typeof,x)
  return_type(f,Ts...)
end

@inline apply_kernel!(::Nothing,f::Function,args...) = f(args...)

# Some particular cases

const NumberOrArray = Union{Number,AbstractArray}

"""
    bcast(f::Function)

Returns a kernel object that represents the "boradcasted" version of the given
function `f`.
"""
bcast(f::Function) = BCasted(f)

struct BCasted{F<:Function} <: Kernel
  f::F
end

function kernel_return_type(f::BCasted,x::Number...)
  Ts = map(typeof,x)
  return_type(f.f,Ts...)
end

function kernel_cache(f::BCasted,x::Number...)
  nothing
end

@inline function apply_kernel!(::Nothing,f::BCasted,x::Number...)
  f.f(x...)
end

function kernel_return_type(f::BCasted,x::NumberOrArray...)
  typeof(kernel_cache(f,x...).array)
end

function kernel_cache(f::BCasted,x::NumberOrArray...)
  s = _sizes(x...)
  bs = Base.Broadcast.broadcast_shape(s...)
  Te = map(numbertype,x)
  T = return_type(f.f,Te...)
  N = length(bs)
  r = Array{T,N}(undef,bs)
  ri = testvalue(T)
  fill!(r,ri)
  cache = CachedArray(r)
   _prepare_cache(cache,x...)
end

numbertype(a::AbstractArray) = eltype(a)

numbertype(a::Number) = typeof(a)

@inline function apply_kernel!(cache,f::BCasted,x::NumberOrArray...)
  r = _prepare_cache(cache,x...)
  a = r.array
  broadcast!(f.f,a,x...)
  a
end

@inline function _prepare_cache(c,x...)
  s = _sizes(x...)
  bs = Base.Broadcast.broadcast_shape(s...)
  if bs != size(c)
    setsize!(c,bs)
  end
  c
end

# TODO use map
@inline function _sizes(a,x...)
  (_size(a), _sizes(x...)...)
end

@inline function _sizes(a)
  (_size(a),)
end

@inline _size(a) = size(a)
@inline _size(a::Number) = (1,)

"""
    elem(f::Function)

Returns a kernel that represents the element-wise
version of the operation `f`
It does not broadcast in singleton axes. Thus, allows some
performance optimizations with respect to broadcast.
"""
elem(f::Function) = Elem(f)

struct Elem{F} <: Kernel
  f::F
  Elem(f::Function) = new{typeof(f)}(f)
end

# TODO more tests here

# It defaults to bcast (TODO test these ones)

@inline function apply_kernel!(cache,k::Elem,x::NumberOrArray...)
  b = bcast(k.f)
  apply_kernel!(cache,b,x...)
end

function kernel_cache(k::Elem,x::NumberOrArray...)
  b = bcast(k.f)
  kernel_cache(b,x...)
end

function kernel_return_type(k::Elem,x::NumberOrArray...)
  b = bcast(k.f)
  kernel_return_type(b,x...)
end

# More Efficient implementations

# Number

function kernel_return_type(k::Elem,a::Number)
  return_type(k.f,typeof(a))
end

function kernel_cache(k::Elem,a::Number)
  nothing
end

@inline function apply_kernel!(c,k::Elem,a::Number)
  k.f(a)
end

# Array

function kernel_return_type(k::Elem,a::AbstractArray)
  typeof(kernel_cache(k,a).array)
end

function kernel_cache(k::Elem,a::AbstractArray)
  T = return_type(k.f,eltype(a))
  CachedArray(similar(a,T))
end

@inline function apply_kernel!(c,f::Elem,a::AbstractArray)
  setsize!(c,size(a))
  r = c.array
  for i in eachindex(a)
    @inbounds r[i] = f.f(a[i])
  end
  r
end

# Array vs Array

function kernel_return_type(k::Elem,a::AbstractArray,b::AbstractArray)
  typeof(kernel_cache(k,a,b).array)
end

function kernel_cache(k::Elem,a::AbstractArray,b::AbstractArray)
  _checks(a,b)
  T = return_type(k.f,eltype(a),eltype(b))
  CachedArray(similar(a,T))
end

@inline function apply_kernel!(c,f::Elem,a::AbstractArray,b::AbstractArray)
  _checks(a,b)
  setsize!(c,size(a))
  r = c.array
  for i in eachindex(a)
    @inbounds r[i] = f.f(a[i],b[i])
  end
  r
end

# Number vs Number

function kernel_return_type(k::Elem,a::Number,b::Number)
  return_type(k.f,typeof(a),typeof(b))
end

function kernel_cache(k::Elem,a::Number,b::Number)
  nothing
end

@inline function apply_kernel!(c,k::Elem,a::Number,b::Number)
  k.f(a,b)
end

# Array vs Number

function kernel_return_type(k::Elem,a::AbstractArray,b::Number)
  typeof(kernel_cache(k,a,b).array)
end

function kernel_cache(k::Elem,a::AbstractArray,b::Number)
  T = return_type(k.f,eltype(a),typeof(b))
  CachedArray(similar(a,T))
end

@inline function apply_kernel!(c,k::Elem,a::AbstractArray,b::Number)
  setsize!(c,size(a))
  r = c.array
  for i in eachindex(a)
    @inbounds r[i] = k.f(a[i],b)
  end
  r
end

# Number vs Array

function kernel_return_type(k::Elem,a::Number,b::AbstractArray)
  typeof(kernel_cache(k,a,b).array)
end

function kernel_cache(k::Elem,a::Number,b::AbstractArray)
  T = return_type(k.f,typeof(a),eltype(b))
  CachedArray(similar(b,T))
end

@inline function apply_kernel!(c,k::Elem,a::Number,b::AbstractArray)
  setsize!(c,size(b))
  r = c.array
  for i in eachindex(b)
    @inbounds r[i] = k.f(a,b[i])
  end
  r
end

function _checks(a,b)
  @assert size(a) == size(b) "Sizes must agree."
  nothing
end

"""
    contract(f::Function)

Like the dot product between to vectors, but using operation `f` instead
of `*` between components.

!!! warning
    not needed any more, to be deleted

# Examples

```jldoctests
using Gridap.Arrays
k = contract(-)
apply_kernel(k,[1,2],[2,4]) # Equivalent to (1-2) + (2-4)
# output
-3
```
"""
contract(f::Function) = Contracted(f)

struct Contracted{F} <: Kernel
  f::F
  Contracted(f::Function) = new{typeof(f)}(f)
end

function kernel_return_type(k::Contracted,a::AbstractArray,b::AbstractArray)
  return_type(k.f,eltype(a),eltype(b))
end

function kernel_cache(k::Contracted,a::AbstractArray,b::AbstractArray)
  kernel_return_type(k,a,b)
end

@inline function apply_kernel!(T,f::Contracted,a::AbstractArray,b::AbstractArray)
  _checks(a,b)
  c = zero(T)
  for i in eachindex(a)
    c += f.f(a[i],b[i])
  end
  c
end

struct MulKernel <: Kernel end

function kernel_cache(k::MulKernel,a,b)
  c = a*b
  CachedArray(c)
end

@inline function apply_kernel!(cache,k::MulKernel,a::AbstractMatrix,b::AbstractVector)
  m = axes(a,1)
  setaxes!(cache,(m,))
  c = cache.array
  mul!(c,a,b)
  c
end

@inline function apply_kernel!(cache,k::MulKernel,a::AbstractMatrix,b::AbstractMatrix)
  m = axes(a,1)
  n = axes(b,2)
  setaxes!(cache,(m,n))
  c = cache.array
  mul!(c,a,b)
  c
end


