"""
Abstract type representing a function (mapping) that provides a cache and an in-place
evaluation for performance. This is the type to be used in the [`apply`](@ref) function.

Derived types must implement the following method:

- [`evaluate!(cache,k,x...)`](@ref)

and optionally these ones:

- [`return_cache(k,x...)`](@ref)
- [`return_type(k,x...)`](@ref)

The mapping interface can be tested with the [`test_mapping`](@ref) function.

Note that most of the functionality implemented in terms of this interface
relies in duck typing. That is, it is not strictly needed to work with types
that inherit from `Mapping`. This is specially useful in order to accommodate
existing types into this framework without the need to implement a wrapper type
that inherits from `Mapping`. For instance, a default implementation is available
for `Function` objects.  However, we recommend that new types inherit from `Mapping`.

"""
abstract type Mapping <: GridapType end

"""
    return_cache(f,x...)

Returns the `cache` needed to apply mapping `f` with arguments
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
evaluate!(cache,f,x...) = @abstractmethod

"""
    return_type(f,x...)

Returns the type of the result of calling mapping `f` with
arguments of the types of the objects `x`.

Its default implementation is `typeof(kernel_testitem(f,x...))`
"""
return_type(f,x...) = typeof(testitem(f,x...))

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

(m::Mapping)(x...) = evaluate(m,x...)

#@fverdugo this default implementation should be improved to be more resilient
# to inputs not defined in the domain of f.
# For instance we could do
#
#     testitem(f,x...) = evaluate(f,testargs(f,x...)...)
#
# with the default implementation
#
#     testargs(f,x...) = x
#
# So that the user can define testargs for functions with non trivial domains:
#
#    myf(x) = sqrt(x-1)
#    testargs(::typeof(myf),x) = zero(x) + 1
#
#
@inline testitem(k,x...) = evaluate(k,x...)

# @fverdugo
# testitem!(cache,k,x...) = evaluate!(cache,k,testargs(k,x...)...)
@inline function testitem!(cache,k,x...)
  evaluate!(cache,k,x...)
end

# Default implementation for Function

evaluate!(cache,f::Function,x...) = f(x...)

# Number or Array implementation
#
# The following default implementation for arrays of number is problematic.
# I would define
#
#    evaluate!(cache,f::GenericFieldArray{AbstractArray{<:Number}},x...) = f.object
#
# instead.
#
# Interpreting arrays as mappings is problematic when you define the apply function
# (see my comments there)
#
# I would try to avoid to interpret any type of array as a mapping.
# In the hypothetical case, we want to give a default mapping behaviour to an array
# I would use this definition: "an array is a mapping from a set of indices to a set of values"
# which is general for any type of array, not only for arrays of numbers. That is:
#
#   return_cache(f::AbstractArray,i...) = array_cache(f,i...)
#   evaluate!(cache,f::AbstractArray,i...) = getindex!(cache,f,i...)
#
# But I believe it is better to not define any default mapping behaviour for arrays of any kind
# since arrays are not callable in Julia.
#
# I also would remove the default Mapping definition for numbers since it is VERY confusing 
#  that evaluate(1,3.0) == 1 and 1(3.0) == 3.0
#  I would define instead
#  evaluate!(cache,a::GenericField{<:Number},x...) = a.object
#

return_type(f::Union{Number,AbstractArray{<:Number}},x...) = typeof(f)

evaluate!(cache,f::Union{Number,AbstractArray{<:Number}},x...) = f


# Testing the interface
"""
    test_mapping(f,x::Tuple,y,cmp=(==))

Function used to test if the mapping `f` has been
implemented correctly. `f` is a `Mapping` sub-type, `x` is a tuple in the domain of the
mapping and `y` is the expected result. Function `cmp` is used to compare
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

# Broadcast Functions

#@fverdugo rename BroadcastMapping -> Broadcasting
"""
    BroadcastMapping(f)

Returns a mapping that represents the "broadcasted" version of the
function `f`.

# Example

```jldoctest
using Gridap.Mappings

a = [3,2]
b = [2,1]

bm = BroadcastMapping(+)

c = evaluate(bm,a,b)

println(c)

# output
[5, 3]
```
"""
struct BroadcastMapping{F} <: Mapping
  f::F
end

# @fverdugo Consider this case:
#     
#     struct Foo end
#     sayhello(a::Foo) = "hi!"
#     @assert sayhello.([Foo(),Foo()]) == ["hi!","hi!"] # Works
#     @assert BroadcastMapping(sayhello)([Foo(),Foo()]) == ["hi!","hi!"]
#
# I would say that previous line does not work with the current implementation.
# We need this more general implementation:
#
#    evaluate!(f::BroadcastMapping,x...) = nothing
#    @inline evaluate!(cache,f::BroadcastMapping,x...) = broadcast(f.f,x...)
#
# and here the last argument should be x::Union{Number,AbstractArray{<:Number}}... also in return_cache below.
@inline function evaluate!(cache,f::BroadcastMapping,x...)
  r = _prepare_cache(cache,x...)
  a = r.array
  broadcast!(f.f,a,x...)
  a
end

function evaluate!(cache,b::BroadcastMapping,args::Number...)
  b.f(args...)
end

function return_type(f::BroadcastMapping,x::Number...)
  Ts = map(typeof,x)
  return_type(f.f,Ts...)
end

function return_type(f::BroadcastMapping,x::AbstractArray...)
  typeof(return_cache(f,x...).array)
end

function return_cache(f::BroadcastMapping,x::Number...)
  nothing
end

#@fverdugo from a so general input, we cannot assume that the result will be an array always
# last argument should be x::Union{Number,AbstractArray{<:Number}}...
function return_cache(f::BroadcastMapping,x...)
  s = _size.(x)
  bs = Base.Broadcast.broadcast_shape(s...)
  Te = map(_numbertype,x)
  T = return_type(f.f,Te...)
  N = length(bs)
  r = Array{T,N}(undef,bs)
  ri = testvalue(T)
  fill!(r,ri)
  cache = CachedArray(r)
  _prepare_cache(cache,x...)
end

@inline _numbertype(a::AbstractArray) = eltype(a)
@inline _numbertype(a::Number) = typeof(a)

@inline function _prepare_cache(c,x...)
  s = _size.(x)
  bs = Base.Broadcast.broadcast_shape(s...)
  if bs != size(c)
    setsize!(c,bs)
  end
  c
end

@inline _size(a) = size(a)
@inline _size(a::Number) = (1,)

# OperationMappings
# @fverdugo I would remove this if it is not used and if we don't expect to use it in the future
# to avoid confusions and to keep focus only in the parts that are used.
"""
    OperationMapping(f,args)

Returns a mapping that represents the result of applying the function `f`
to the arguments in the tuple `args`.
"""
struct OperationMapping{K,L} <: Mapping
  k::K
  l::L
  @inline function OperationMapping(k,l)
    new{typeof(k),typeof(l)}(k,l)
  end
end

function return_type(c::OperationMapping,x...)
  Ts = map(fi -> return_type(fi,x...),c.l)
  return_type(c.k, testvalues(Ts...)...)
end

function return_cache(c::OperationMapping,x...)
  cl = map(fi -> return_cache(fi,x...),c.l)
  lx = map((ci,fi) -> evaluate!(ci,fi,x...),cl,c.l)
  ck = return_cache(c.k,lx...)
  (ck,cl)
end

@inline function evaluate!(cache,c::OperationMapping,x...)
  ck, cf = cache
  lx = map((ci,fi) -> evaluate!(ci,fi,x...),cf,c.l)
  evaluate!(ck,c.k,lx...)
end

# Operations

#@fverdugo The result is not always an OperationMapping, so I would not state this in the documentation.
# In practice it will be an OperationField.
# I would remove OperationMapping and move Operation to the source file
# where operation between fields are defined (and document only
# Operation for Field arguments). I think this would be much more clear to understand for new users.
"""
    Operation(op)

Returns a mapping that, when applied to a tuple `args`, returns
`OperationMapping(op,args)`

# Example

```jldoctest
using Gridap.Mappings

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
struct Operation{T} <: Mapping
  op::T
end

"""
    operation(op)

    Idem as [`Operation(op)`](@ref)).
"""
operation(a) = Operation(a)

evaluate!(cache,op::Operation,x...) = OperationMapping(op.op,x)

# @fverdugo
# this is not needed since Operation <: Mapping and functor-like evaluation
# is already defined for Mapping
(op::Operation)(x...) = evaluate!(nothing,op,x...)
