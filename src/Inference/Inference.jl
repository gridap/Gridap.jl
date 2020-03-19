"""
This module provides a set of helper function to safely infer return types of functions.

In Gridap, we rely as less as possible in type inference. But, when needed, we adopt
the following mechanism in order to compute returned types. We do not rely on
the `Base._return_type` function.

This module exports following functions:
$(EXPORTS)

"""
module Inference

using DocStringExtensions
using FillArrays
using Gridap.Helpers

export testvalue
export testvalues
export testargs
export testargs_broadcast
export return_type
export return_type_broadcast

"""
$(TYPEDSIGNATURES)

Returns the type returned by function `f` when called with arguments
of the types in `Ts`.

The underlying implementation uses the function [`testargs`](@ref) to generate
some test values in order to call the function and determine the returned type.
This mechanism does not use `Base._return_type`. One of the advantages is
that the given function `f` is called, and thus, meaningful error messages
will be displayed if there is any error in `f`.

"""
function return_type(f::Function,Ts...)
  args = testargs(f,Ts...)
  #try
    typeof(f(args...))
  #catch e
  #  if isa(e,DomainError)
  #    s = "Function $(nameof(f)) cannot be evaluated at $args, its not in the domain.\n"
  #    s *= " Define function `testargs(::typeof{$(nameof(f))},Ts...)`\n"
  #    s *= " which sould return an argument tuple in the function domain."
  #    error(s)
  #  else
  #    throw(e)
  #  end
  #end
  # TODO try - catch block makes function very in-efficient when
  # called multiple times
end

"""
    return_type_broadcast(f::Function,Ts::DataType...) -> DataType

Like [`return_type`](@ref), but when function `f` is used in a broadcast operation.
"""
function return_type_broadcast(f::Function,Ts...)
  args = testargs_broadcast(f,Ts...)
  r = broadcast(f,args...)
  typeof(r)
end

function return_type_broadcast(f::Function,Ts::Type{<:Number}...)
  return_type(f,Ts...)
end

"""
$(SIGNATURES)
"""
function testargs_broadcast(f::Function,Ts...)
  testvalues(Ts...)
end

function testargs_broadcast(f::Function,Ts::Type{<:Number}...)
  testargs(f,Ts...)
end

"""
    testargs(f::Function,Ts::DataType...) -> Tuple

Returns a tuple with valid arguments of the types in `Ts` in order to call
function `f`. It defaults to `testvalues(Ts...)`, see the [`testvalues`](@ref)
function.
The user can overload the `testargs`
function for particular functions if the default test arguments are not in the domain
of the function and a `DomainError` is raised.

# Examples

For the following function, the default test argument (which is a zero)
is not in the domain. We can overload the `testargs` function to provide
a valid test argument.

```jldoctests
using Gridap.Inference
import Gridap.Inference: testargs
foo(x) = sqrt(x-1)
testargs(::typeof(foo),T::DataType) = (one(T),)
return_type(foo, Int)
# output
Float64
```

"""
testargs(f::Function,Ts...) = testvalues(Ts...)

"""
    testvalue(::Type{T}) where T

Returns an arbitrary instance of type `T`. It defaults to `zero(T)` for
non-array types and to an empty array for array types.
This function is used to compute the default test arguments in
[`testargs`](@ref).
It can be overloaded for new types `T` if `zero(T)` does not makes sense.
"""
function testvalue end

testvalue(::Type{T}) where T = zero(T)
testvalue(v) = testvalue(typeof(v))

function testvalue(::Type{T}) where T<:AbstractArray{E,N} where {E,N}
   similar(T,fill(0,N)...)
end

function testvalue(::Type{T}) where T<:Fill{E,N} where {E,N}
  Fill(zero(E),fill(0,N)...)
end

function testvalue(::Type{<:Tuple})
  @notimplemented "testvalue on Tuple type only implemented up to 8 tuple elements"
end

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
    testvalues(Ts::DataType...) -> Tuple

Returns a tuple with test values for each of the types in `Ts`.
Equivalent to `map(testvalue,Ts)`.
"""
function testvalues(a,b...)
  ta = testvalue(a)
  tb = testvalues(b...)
  (ta,tb...)
end

function testvalues(a)
  ta = testvalue(a)
  (ta,)
end

end # module
