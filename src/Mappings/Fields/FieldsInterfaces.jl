"""
    const Point{D,T} = VectorValue{D,T}

Type representing a point of D dimensions with coordinates of type T.
Fields are evaluated at vectors of `Point` objects.
"""
const Point{D,T} = VectorValue{D,T}

"""
    abstract type Field <: Mapping

Abstract type representing physical fields, bases of fields, and other related objects.
These different cases are distinguished by the return value obtained when evaluating them. E.g.,
a physical field returns a vector of values when evaluated at a vector of points, and a basis of `nf` fields
returns a 2d matrix (`np` x `nf`) when evaluated at a vector of `np` points.

The following functions need to be overloaded:

- [`evaluate_field!(cache,f,x)`](@ref)
- [`field_cache(f,x)`](@ref)

The following functions can be also provided optionally

- [`field_gradient(f)`](@ref)
- [`field_return_type(f,x)`](@ref)

Moreover, if the [`field_gradient(f)`](@ref) is not provided, a default implementation that uses the following
functions will be used.

- [`evaluate_gradient!(cache,f,x)`](@ref)
- [`gradient_cache(f,x)`](@ref)

In order to be able to call `field_gradient` again on the resulting object the following methods have to be
provided

- [`evaluate_hessian!(cache,f,x)`](@ref)
- [`hessian_cache(f,x)`](@ref)

These four methods are only designed to be called by the default implementation of [`field_gradient(f)`](@ref) and thus
cannot be assumed that they are available for an arbitrary field. For this reason, these functions are not
exported. The general way of evaluating a gradient of a field is to
build the gradient with [`field_gradient(f)`](@ref) and evaluating the resulting object. For evaluating
the hessian, use two times `field_gradient`.

The interface can be tested with

- [`test_field`](@ref)

Most of the functionality implemented in terms of this interface relies in duck typing (this is why all functions in the interface
have the word "field").  Thus, it is not strictly needed to work with types
that inherit from `Field`. This is specially useful in order to accommodate
existing types into this framework without the need to implement a wrapper type that inherits from `Field`.
For instance, a default implementation is available for numbers, which behave like "constant" fields, or arrays of numbers, which behave like
"constant" bases of fields.  However, we recommend that new types inherit from `Field`.

"""
abstract type NewField <: Mapping end

const FieldVector{T} = AbstractVector{T} where T<:NewField

const FieldArray{T,N} = AbstractArray{T,N} where {T<:NewField,N}

const FieldOrFieldArray = Union{NewField,FieldArray}

function return_type(f::FieldOrFieldArray,x::AbstractArray{<:Point})
  typeof(evaluate(f,x))
end

# Default implementation for the non-vectorised case

return_type(f::NewField,x::Point) = return_type(f,Fill(x,1))

function evaluate!(cache,f::NewField,x::Point)
  first(evaluate!(cache,f,Fill(x,1)))
end

function evaluate(f::NewField,x::Point)
  first(evaluate(f,Fill(x,1)))
end

function evaluate!(cache,f::AbstractArray{<:NewField},x::Point)
  @notimplemented
end

"""
    gradient_type(::Type{T},x::Point) where T
"""
function return_gradient_type(::Type{T},x::Point) where T
  typeof(outer(zero(x),zero(T)))
end

function gradient end

@inline derivative(f::NewField) = transpose(gradient(f))

"""
    const ∇ = gradient

Alias for the `gradient` function.
"""

function evaluate_gradient!(cache,f,x)
  @abstractmethod
end

# """
# # $(SIGNATURES)
# # """
function return_gradient_cache(f,x)
#   @abstractmethod
end

# """
# $(SIGNATURES)
# """
function evaluate_hessian!(cache,f,x)
  @abstractmethod
end

# """
# $(SIGNATURES)
# """
function return_hessian_cache(f,x)
  @abstractmethod
end
# Testers

"""
    test_field(
      f,
      x::AbstractVector{<:Point},
      v::AbstractArray,cmp=(==);
      grad=nothing,
      hessian=nothing)

Function used to test the field interface. `v` is an array containing the expected
result of evaluating the field `f` at the vector of points `x`. The comparison is performed using
the `cmp` function. For fields objects that support the `field_gradient` function, the key-word
argument `grad` can be used. It should contain the result of evaluating `field_gradient(f)` at x.
Idem for `hessian`.
The checks are performed with the `@test` macro.
"""
function test_field(
  f::FieldOrFieldArray,
  x::Tuple,
  v::AbstractArray,cmp=(==);
  grad=nothing,
  hessian=nothing)

  x, = x

  @test isa(x,AbstractVector{<:Point})

  w = evaluate(f,x)

  np, = size(w)
  @test length(x) == np
  @test cmp(w,v)
  @test typeof(w) == return_type(f,x)

  cf = return_cache(f,x)
  r = evaluate!(cf,f,x)
  @test cmp(r,v)

  _x = vcat(x,x)
  _v = vcat(v,v)
  _w = evaluate!(cf,f,_x)
  @test cmp(_w,_v)

  if isa(f,NewField)
    test_mapping(f,(x,),v,cmp)
  end

  if grad != nothing
    g = gradient(f)
    if typeof(f) <: NewField
      @test g isa NewField
    elseif typeof(f) <: AbstractArray{<:NewField}
      @test g isa AbstractArray{<:NewField}
    end
    test_field(g,(x,),grad,cmp,grad=hessian)
  end

end

function test_field(f::FieldOrFieldArray,x,v::AbstractArray,cmp=(==);grad=nothing,hessian=nothing)
  test_field(f,(x,),v,cmp;grad=grad,hessian=hessian)
end
