
"""
    abstract type Dof <: Mapping

Abstract type representing a degree of freedom (DOF), a basis of DOFs, and related objects.
These different cases are distinguished by the return type obtained when evaluating the `Dof`
object on a `Field` object. See function [`evaluate_dof!`](@ref) for more details.

The following functions needs to be overloaded

- [`dof_cache`](@ref)
- [`evaluate_dof!`](@ref)

The following functions can be overloaded optionally

- [`dof_return_type`](@ref)

The interface is tested with

- [`test_dof`](@ref)

In most of the cases it is not strictly needed that types that implement this interface
inherit from `Dof`. However, we recommend to inherit from `Dof`, when possible.


"""
abstract type Dof <: Mapping end

"""
    dof_cache(dof,field)

Returns the cache needed to call `evaluate_dof!(cache,dof,field)`
"""
function dof_cache(dof,field)
  @abstractmethod
end

"""
    evaluate_dof!(cache,dof,field)

Evaluates the dof `dof` with the field `field`. It can return either an scalar value or
an array of scalar values depending the case. The `cache` object is computed with function
[`dof_cache`](@ref).

When a mathematical dof is evaluated on a physical field, a scalar number is returned. If either
the `Dof` object is a basis of DOFs, or the `Field` object is a basis of fields,
or both objects are bases, then the returned object is an array of scalar numbers. The first
dimensions in the resulting array are for the `Dof` object and the last ones for the `Field`
object. E.g, a basis of `nd` DOFs evaluated at physical field returns a vector of `nd` entries.
A basis of `nd` DOFs evaluated at a basis of `nf` fields returns a matrix of size `(nd,nf)`.
"""
function evaluate_dof!(cache,dof,field)
  @abstractmethod
end

"""
    dof_return_type(dof,field)

Returns the type for the value obtained with evaluating `dof` with `field`.

It defaults to

    typeof(evaluate_dof(dof,field))
"""
function dof_return_type(dof,field)
  typeof(evaluate_dof(dof,field))
end

# Testers

"""
    test_dof(dof,field,v,comp::Function=(==))

Test that the `Dof` interface is properly implemented
for object `dof`. It also checks if the object `dof`
when evaluated at the field `field` returns the same
value as `v`. Comparison is made with the `comp` function.
"""
function test_dof(dof,field,v,comp::Function=(==))
  if isa(dof,Dof)
    test_mapping(dof,(field,),v,comp)
  end
  r = evaluate_dof(dof,field)
  @test comp(r,v)
  @test typeof(r) == dof_return_type(dof,field)
end

# Implement Mapping interface

@inline return_cache(dof::Dof,field) = dof_cache(dof,field)

@inline evaluate!(cache,dof::Dof,field) = evaluate_dof!(cache,dof,field)

@inline return_type(dof::Dof,field) = dof_return_type(dof,field)

# Some API

"""
    evaluate_dof(dof,field)

Equivalent to

    cache = dof_cache(dof,field)
    evaluate_dof!(cache,dof,field)
"""
function evaluate_dof(dof,field)
  cache = dof_cache(dof,field)
  evaluate_dof!(cache,dof,field)
end

"""
    evaluate(dof::Dof,field)

Equivalent to `evaluate_dof(dof,field)`.
"""
evaluate(dof::Dof,field) = evaluate_dof(dof,field)

# Working with arrays of Dofs

"""
    evaluate_dof_array(dof::AbstractArray,field::AbstractArray)

Evaluates the `Dof` objects in the array `dof` at the `Field` objects
at the array `field` element by element.

The result is numerically equivalent to

    map(evaluate_dof, dof, field)

but it is described with a more memory-friendly lazy type.
"""
function evaluate_dof_array(dof::AbstractArray,field::AbstractArray)
  k = DofEval()
  lazy_map(k,dof,field)
end

function evaluate_dof_array(dof::AbstractArray{<:Dof},field::AbstractArray)
  lazy_map(dof,field)
end

"""
    evaluate(dof::AbstractArray{<:Dof},field::AbstractArray)

Equivalent to `evaluate_dof_array(dof,field)`
"""
function evaluate(dof::AbstractArray{<:Dof},field::AbstractArray)
  evaluate_dof_array(dof,field)
end

struct DofEval <: Mapping end

function return_cache(k::DofEval,dof,field)
  dof_cache(dof,field)
end

@inline function evaluate!(cache,k::DofEval,dof,field)
  evaluate_dof!(cache,dof,field)
end

function return_type(k::DofEval,dof,field)
  dof_return_type(dof,field)
end
