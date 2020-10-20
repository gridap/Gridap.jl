"""
    const Point{D,T} = VectorValue{D,T}

Type representing a point of D dimensions with coordinates of type T.
Fields are evaluated at vectors of `Point` objects.
"""
const Point{D,T} = VectorValue{D,T}

"""
    abstract type Field <: Map

Abstract type representing a physical (scalar, vector, or tensor) field. The
domain is a `Point` and the range a scalar (i.e., a sub-type of Julia `Number`),
a `VectorValue`, or a `TensorValue`.

These different cases are distinguished by the return value obtained when evaluating them. E.g.,
a physical field returns a vector of values when evaluated at a vector of points, and a basis of `nf` fields
returns a 2d matrix (`np` x `nf`) when evaluated at a vector of `np` points.

The following functions (i.e., the `Map` API) need to be overloaded:

- [`evaluate!(cache,f,x)`](@ref)
- [`return_cache(f,x)`](@ref)

and optionally

- [`return_type(f,x)`](@ref)

A `Field` can also provide its gradient if the following function is implemented
- [`gradient(f)`](@ref)

Higher derivatives can be obtained if the resulting object also implements this method.

Moreover, if the [`gradient(f)`](@ref) is not provided, a default implementation that uses the
following functions will be used.

- [`evaluate_gradient!(cache,f,x)`](@ref)
- [`return_gradient_cache(f,x)`](@ref)

Higher order derivatives require the implementation of

- [`evaluate_hessian!(cache,f,x)`](@ref)
- [`return_hessian_cache(f,x)`](@ref)

These four methods are only designed to be called by the default implementation of [`field_gradient(f)`](@ref) and thus
cannot be assumed that they are available for an arbitrary field. For this reason, these functions are not
exported. The general way of evaluating a gradient of a field is to
build the gradient with [`gradient(f)`](@ref) and evaluating the resulting object. For evaluating
the hessian, use two times `gradient`.

The interface can be tested with

- [`test_field`](@ref)

For performance, the user can also consider a _vectorised_ version of the
`Field` API that evaluates the field in a vector of points (instead of only one
point). E.g., the `evaluate!` function for a vector of points returns a vector
of scalar, vector or tensor values.

"""
abstract type Field <: Map end

evaluate!(c,f::Field,x::Point) = @abstractmethod

# Differentiation

function gradient end
hessian(f) = gradient(gradient(f))

const ∇ = gradient
const ∇∇ = hessian


# Hook methods used by FieldGradient and FieldHessian

return_gradient_cache(f::Field,x::Point) = nothing

evaluate_gradient!(c,f::Field,x::Point) = @abstractmethod

testargs_gradient(f::Field,x::Point) = testargs(f,x)

function return_gradient_value(f::Field,x::Point)
  y = testargs_gradient(f,x)
  c = return_gradient_cache(f,y...)
  evaluate_gradient!(c,f,y...)
end

return_hessian_cache(f::Field,x::Point) = nothing

evaluate_hessian!(c,f::Field,x::Point) = @abstractmethod

testargs_hessian(f::Field,x::Point) = testargs(f,x)

function return_hessian_value(f::Field,x::Point)
  y = testargs_hessian(f,x)
  c = return_hessian_cache(f,y...)
  evaluate_hessian!(c,f,y...)
end

# Default methods for arrays of points

function testargs(f::Field,x::AbstractArray{<:Point})
  y = copy(x)
  broadcast!(xi->first(testargs(f,xi)),y,x)
  (y,)
end

@inline function return_cache(f::Field,x::AbstractArray{<:Point})
  T = return_type(f,testitem(x))
  s = size(x)
  ab = zeros(T,s)
  cb = CachedArray(ab)
  cf = return_cache(f,testitem(x))
  cb, cf
end

@inline function evaluate!(c,f::Field,x::AbstractArray{<:Point})
  cb, cf = c
  sx = size(x)
  setsize!(cb,sx)
  r = cb.array
  for i in eachindex(x)
    @inbounds r[i] = evaluate!(cf,f,x[i])
  end
  r
end

# Default methods for arrays of points (gradient)

@inline function return_gradient_cache(f::Field,x::AbstractArray{<:Point})
  cf = return_gradient_cache(f,testitem(x))
  T = typeof(return_gradient_value(f,testitem(x)))
  s = size(x)
  ab = zeros(T,s)
  cb = CachedArray(ab)
  cb, cf
end

@inline function evaluate_gradient!(c,f::Field,x::AbstractArray{<:Point})
  cb, cf = c
  sx = size(x)
  setsize!(cb,sx)
  r = cb.array
  for i in eachindex(x)
    @inbounds r[i] = evaluate_gradient!(cf,f,x[i])
  end
  r
end

# Default methods for arrays of points (hessian)

@inline function return_hessian_cache(f::Field,x::AbstractArray{<:Point})
  cf = return_hessian_cache(f,testitem(x))
  T = typeof(return_hessian_value(f,testitem(x)))
  s = size(x)
  ab = zeros(T,s)
  cb = CachedArray(ab)
  cb, cf
end

@inline function evaluate_hessian!(c,f::Field,x::AbstractArray{<:Point})
  cb, cf = c
  sx = size(x)
  setsize!(cb,sx)
  r = cb.array
  for i in eachindex(x)
    @inbounds r[i] = evaluate_hessian!(cf,f,x[i])
  end
  r
end

# GenericField

"""
A wrapper for objects that can act as fields, e.g., functions which implement the `Field` API.
"""
struct GenericField{T} <: Field
  object::T
end

#@inline Field(f) = GenericField(f)
@inline GenericField(f::Field) = f

testargs(a::GenericField,x::Point) = testargs(a.object,x)
return_value(a::GenericField,x::Point) = return_value(a.object,x)
return_cache(a::GenericField,x::Point) = return_cache(a.object,x)
@inline evaluate!(cache,a::GenericField,x::Point) = evaluate!(cache,a.object,x)
gradient(a::GenericField) = GenericField(gradient(a.object))

#@inline return_type(::Type{<:GenericField},::T) where T<:Field = T
#@inline return_type(::Type{<:GenericField},::T) where T = GenericField{T}
#@inline return_type(a::GenericField,x) = return_type(a.object,x)

# Make Field behave like a collection

@inline Base.length(::Field) = 1
@inline Base.size(::Field) = ()
@inline Base.axes(::Field) = ()
@inline Base.IteratorSize(::Type{<:Field}) = Base.HasShape{0}()
@inline Base.eltype(::Type{T}) where T<:Field = T
@inline Base.iterate(a::Field) = (a,nothing)
@inline Base.iterate(a::Field,::Nothing) = nothing
@inline Base.getindex(a::Field,i::Integer) =  (@check i == 1; a)

# Zero field

@inline Base.zero(a::Field) = ZeroField(a)

"""
It represents `0.0*f` for a field `f`.
"""
struct ZeroField{F} <: Field
  field::F
end

return_cache(z::ZeroField,x::Point) = zero(return_type(z.field,x))
@inline evaluate!(cache,z::ZeroField,x::Point) = cache

function return_cache(z::ZeroField,x::AbstractArray{<:Point})
  E = return_type(z.field,testitem(x))
  c = zeros(E,size(x))
  CachedArray(c)
end

function evaluate!(c,f::ZeroField,x::AbstractArray{<:Point})
  nx = size(x)
  if size(c) != nx
    setsize!(c,nx)
    fill!(c.array,zero(eltype(c)))
  end
  c.array
end

# @inline function evaluate_gradient!(cache,z::ZeroField,x::Point)
#   outer(zero(return_type(z.field)),zero(x))
# end

# @inline function evaluate_gradient!(cache,z::ZeroField,x::AbstractArray{<:Point})
#   outer(zero(return_type(z.field)),zero(x))
# end

@inline gradient(z::ZeroField) = ZeroField(gradient(z.field))

# Make Number behave like Field
#

#const ConstantField{T} = GenericField{T} where T<:Number

# Number itself does not implement the Field interface since Number objects are not callable in Julia.
# Thus, wrapping a Number in a GenericField does not makes sense since the wrapped object
# is assumed to implement the Field interface.
# I think it is conceptually better to have ConstantField as struct otherwise we break the invariant
# "for any object wrapped in a GenericField we can assume that it implements the Field interface"
struct ConstantField{T<:Number} <: Field
  object::T
end

#@inline testargs(a::ConstantField,x) = (x,)
#@inline return_value(a::ConstantField,x) = evaluate(a,x)
#@inline return_type(a::ConstantField,x) = typeof(return_value(a,x))

#@inline return_type(f::ConstantField,x) = typeof(f.object)

@inline function evaluate!(c,f::ConstantField,x::Point)
  f.object
end

#function return_type(f::ConstantField,x::AbstractArray{<:Point})
#  typeof(return_cache(f,x).array)
#end

function return_cache(f::ConstantField,x::AbstractArray{<:Point})
  nx = size(x)
  c = fill(f.object,nx)
  CachedArray(c)
end

function evaluate!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = size(x)
  if size(c) != nx
    setsize!(c,nx)
    fill!(c.array,f.object)
  end
  c.array
end

@inline function return_gradient_cache(f::ConstantField,x::Point)
  gradient(f.object)(x)
end

@inline evaluate_gradient!(c,f::ConstantField,x::Point) = c

@inline function return_gradient_cache(f::ConstantField,x::AbstractArray{<:Point})
  CachedArray(gradient(f.object).(x))
end

function evaluate_gradient!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = size(x)
  if size(c) != nx
    setsize!(c,nx)
    fill!(c.array,zero(eltype(c)))
  end
  c.array
end

@inline function return_hessian_cache(f::ConstantField,x::Point)
  hessian(f.object)(x)
end

@inline function return_hessian_cache(f::ConstantField,x::AbstractArray{<:Point})
  CachedArray(hessian(f.object).(x))
end

@inline evaluate_hessian!(c,f::ConstantField,x::Point) = c

function evaluate_hessian!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = size(x)
  if size(c) != nx
    setsize!(c,nx)
    fill!(c.array,zero(eltype(c)))
  end
  c.array
end

# Not needed any more
## Make Function behave like Field
#
#const FunctionField{F} = GenericField{F} where F<:Function
#
##function return_cache(f::FunctionField,x::AbstractArray{<:Point})
#  nx = length(x)
#  c = zeros(return_type(f.object,testitem(x)),nx)
#  CachedArray(c)
#end
#
#function return_type(f::FunctionField,x::AbstractArray{<:Point})
#  typeof(return_cache(f,x).array)
#end
#
#function evaluate!(c,f::FunctionField,x::AbstractArray{<:Point})
#  nx = length(x)
#  setsize!(c,(nx,))
#  for i in eachindex(x)
#    c[i] = f.object(x[i])
#  end
#  c.array
#end
#
#@inline function return_gradient_cache(f::FunctionField,x::Point)
#  gradient(f.object)
#
#end
#
#function return_gradient_cache(f::FunctionField,x::AbstractArray{<:Point})
#  gf = gradient(f.object)
#  nx = length(x)
#  c = zeros(return_type(gf,testitem(x)),nx)
#  # gf, CachedArray(c)
#  gf, CachedArray(c)
#end
#
#@inline evaluate_gradient!(c,f::FunctionField,x::Point) = c(x)
#
#function evaluate_gradient!(cache,f::FunctionField,x::AbstractArray{<:Point})
#  gf, c = cache
#  nx = length(x)
#  setsize!(c,(nx,))
#  for i in eachindex(x)
#    c[i] = gf(x[i])
#  end
#  c.array
#end

# Differentiation

"""
Type that represents the gradient of a field. The wrapped field implements must
implement `evaluate_gradient!` and `return_gradient_cache` for this gradient
to work.
"""
struct FieldGradient{F} <: Field
  object::F
end

gradient(f::Field) = FieldGradient(f)

#@inline gradient(f::GenericField{FieldGradient}) = FieldHessian(f.object.object)

testargs(f::FieldGradient,x::Point) = testargs_gradient(f.object,x)
return_value(f::FieldGradient,x::Point) = return_gradient_value(f.object,x)
return_cache(f::FieldGradient,x::Point) = return_gradient_cache(f.object,x)
@inline evaluate!(cache,f::FieldGradient,x::Point) = evaluate_gradient!(cache,f.object,x)

return_cache(f::FieldGradient,x::AbstractArray{<:Point}) = return_gradient_cache(f.object,x)
@inline evaluate!(cache,f::FieldGradient,x::AbstractArray{<:Point}) = evaluate_gradient!(cache,f.object,x)

"""
Type that represents the hessian of a field. The wrapped field implements must
implement `evaluate_hessian!` and `return_hessian_cache` for this Hessian
to work.
"""
struct FieldHessian{F} <: Field
  object::F
end

gradient(f::FieldGradient) = FieldHessian(f.object)

@inline function gradient(f::FieldHessian)
  @unreachable "Default implementation of 3rt order derivatives not available"
end

testargs(f::FieldHessian,x::Point) = testargs_hessian(f.object,x)
return_value(f::FieldHessian,x::Point) = return_hessian_value(f.object,x)
return_cache(f::FieldHessian,x::Point) = return_hessian_cache(f.object,x)
@inline evaluate!(cache,f::FieldHessian,x::Point) = evaluate_hessian!(cache,f.object,x)

return_cache(f::FieldHessian,x::AbstractArray{<:Point}) = return_hessian_cache(f.object,x)
@inline evaluate!(cache,f::FieldHessian,x::AbstractArray{<:Point}) = evaluate_hessian!(cache,f.object,x)

# Operations

"""
A `Field` that is obtained as a given operation over a tuple of fields.
"""
struct OperationField{O,F} <: Field
  op::O
  fields::F
end

function return_value(c::OperationField,x::Point)
  fx = map(f -> return_value(f,x),c.fields)
  return_value(c.op,fx...)
end

function return_cache(c::OperationField,x::Point)
  cl = map(fi -> return_cache(fi,x),c.fields)
  lx = map(fi -> return_value(fi,x),c.fields)
  ck = return_cache(c.op,lx...)
  ck, cl
end

@inline function evaluate!(cache,c::OperationField,x::Point)
  ck, cf = cache
  lx = map((ci,fi) -> evaluate!(ci,fi,x),cf,c.fields)
  evaluate!(ck,c.op,lx...)
end

function return_cache(c::OperationField,x::AbstractArray{<:Point})
  cf = map(fi -> return_cache(fi,x),c.fields)
  lx = map((ci,fi) -> evaluate!(ci,fi,x),cf,c.fields)
  ck = return_cache(c.op,map(testitem,lx)...)
  r = c.op.(lx...)
  ca = CachedArray(r)
  ca, ck, cf
end

@inline function evaluate!(cache,c::OperationField,x::AbstractArray{<:Point})
  ca, ck, cf = cache
  sx = size(x)
  setsize!(ca,sx)
  lx = map((ci,fi) -> evaluate!(ci,fi,x),cf,c.fields)
  r = ca.array
  for i in eachindex(x)
    @inbounds r[i] = evaluate!(ck,c.op,map(lxi -> lxi[i], lx)...)
  end
  r
end

@inline evaluate!(cache,op::Operation,x::Field...) = OperationField(op.op,x)

for op in (:+,:-,:*,:⋅,:inv,:det)
  @eval ($op)(a::Field...) = Operation($op)(a...)
end

# Operation rules e.g.
for op in (:+,:-)
  @eval begin
    function gradient(a::OperationField{typeof($op)})
      f = a.fields
      g = map( gradient, f)
      $op(g...)
    end
  end
end

# Some syntactic sugar
@inline *(A::Number, B::Field) = ConstantField(A)*B
@inline *(A::Field, B::Number) = ConstantField(B)*A
@inline *(A::Function, B::Field) = GenericField(A)*B
@inline *(A::Field, B::Function) = GenericField(B)*A

# Product rules for the gradient

function product_rule(fun,f1,f2,∇f1,∇f2)
  msg = "Product rule not implemented for product $fun between types $(typeof(f1)) and $(typeof(f2))"
  @notimplemented msg
end

function product_rule(::typeof(*),f1::Real,f2::Real,∇f1,∇f2)
  ∇f1*f2 + f1*∇f2
end

function product_rule(::typeof(*),f1::Real,f2::VectorValue,∇f1,∇f2)
  ∇f1⊗f2 + ∇f2*f1
end

function product_rule(::typeof(*),f1::VectorValue,f2::Real,∇f1,∇f2)
  product_rule(*,f2,f1,∇f2,∇f1)
end

function product_rule(::typeof(⋅),f1::VectorValue,f2::VectorValue,∇f1,∇f2)
  ∇f1⋅f2 + ∇f2⋅f1
end

for op in (:*,:⋅,:⊙,:⊗)
  @eval begin
    function gradient(a::OperationField{typeof($op)})
      f = a.fields
      @notimplementedif length(f) != 2
      f1, f2 = f
      g1, g2 = map(gradient, f)
      k(F1,F2,G1,G2) = product_rule($op,F1,F2,G1,G2)
      Operation(k)(f1,f2,g1,g2)
    end
  end
end

# Chain rule
function gradient(f::OperationField{<:Field})
  a = f.op
  @notimplementedif length(f.fields) != 1
  b, = f.fields
  x = ∇(a)∘b
  y = ∇(b)
  y⋅x
end

# Composition

"""
    f∘g

It returns the composition of two fields, which is just `Operation(f)(g)`
"""
@inline Base.:∘(f::Field,g::Field) = Operation(f)(g)

# Other operations

@inline transpose(f::Field) = f
#@inline Base.copy(f::Field) = f
@inline *(f::Field,g::Field) = Operation(*)(f,g)#⋅g

# Testers

"""
    test_field(
      f::Union{Field,AbstractArray{<:Field}},
      x,
      v,
      cmp=(==);
      grad=nothing,
      hessian=nothing)

Function used to test the field interface. `v` is an array containing the expected
result of evaluating the field `f` at the point or vector of points `x`. The comparison is performed using
the `cmp` function. For fields objects that support the `gradient` function, the keyword
argument `grad` can be used. It should contain the result of evaluating `gradient(f)` at x.
Idem for `hessian`. The checks are performed with the `@test` macro.
"""
function test_field(f::Field, x, v, cmp=(==); grad=nothing, hess=nothing)
  test_mapping(f,(x,),v,cmp)
  if grad != nothing
    test_mapping(gradient(f),(x,),grad,cmp)
  end
  if hess != nothing
    test_mapping(hessian(f),(x,),hess,cmp)
  end
end




#function test_field(
#  f::Union{Field,AbstractArray{<:Field}},
#  x::Tuple,
#  v,
#  cmp=(==);
#  grad=nothing,
#  hessian=nothing)
#
#  x, = x
#
#  @test isa(x,Union{Point,AbstractArray{<:Point}})
#
#  w = evaluate(f,x)
#
#  @test cmp(w,v)
#  @test typeof(w) == return_type(f,x)
#
#  cf = return_cache(f,x)
#  r = evaluate!(cf,f,x)
#  @test cmp(r,v)
#
#  if x isa AbstractArray{<:Point}
#
#    _x = vcat(x,x)
#    _v = vcat(v,v)
#    _w = evaluate!(cf,f,_x)
#    @test cmp(_w,_v)
#  end
#
#  if isa(f,Field)
#    test_mapping(f,(x,),v,cmp)
#  end
#
#  if grad != nothing
#    g = gradient(f)
#    if typeof(f) <: Field
#      @test g isa Field
#    elseif typeof(f) <: AbstractArray{<:Field}
#      @test g isa AbstractArray{<:Field}
#    end
#    test_field(g,(x,),grad,cmp,grad=hessian)
#  end
#end
#
#@inline function test_field(f::Union{Field,AbstractArray{<:Field}},x,v,cmp=(==);grad=nothing,hessian=nothing)
#  test_field(f,(x,),v,cmp;grad=grad,hessian=hessian)
#end
