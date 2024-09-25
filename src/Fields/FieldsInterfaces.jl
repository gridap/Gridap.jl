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


The next paragraph is out-of-date:

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
const ∇ = gradient
∇∇(f) = gradient(gradient(f))

gradient(f,::Val{1}) = ∇(f)
gradient(f,::Val{2}) = ∇∇(f)

evaluate!(cache,::Broadcasting{typeof(∇)},a::Field) = ∇(a)
evaluate!(cache,::Broadcasting{typeof(∇∇)},a::Field) = ∇∇(a)
lazy_map(::Broadcasting{typeof(∇)},a::AbstractArray{<:Field}) = lazy_map(∇,a)
lazy_map(::Broadcasting{typeof(∇∇)},a::AbstractArray{<:Field}) = lazy_map(∇∇,a)

push_∇(∇a::Field,ϕ::Field) = pinvJt(∇(ϕ))⋅∇a

function pinvJt(Jt::MultiValue{Tuple{D,D}}) where D
  inv(Jt)
end

function pinvJt(Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  @check D1 < D2
  J = transpose(Jt)
  transpose(inv(Jt⋅J)⋅Jt)
end

function push_∇∇(∇∇a::Field,ϕ::Field)
  @notimplemented """\n
  Second order derivatives of quantities defined in the reference domain not implemented yet.

  This is a feature that we want to have at some point in Gridap.
  If you are ready to help with this implementation, please contact the
  Gridap administrators.
  """
end

"""
    gradient_type(::Type{T},x::Point) where T
"""
function gradient_type(::Type{T},x::Point) where T
  typeof(outer(zero(x),zero(T)))
end

"""
Type that represents the gradient of a field. The wrapped field must
implement `evaluate_gradient!` and `return_gradient_cache` for this gradient
to work.

N is how many times the gradient is applied
"""
struct FieldGradient{N,F} <: Field
  object::F
  FieldGradient{N}(object::F) where {N,F} = new{N,F}(object)
end

gradient(f::Field) = FieldGradient{1}(f)
gradient(f::FieldGradient{N}) where N = FieldGradient{N+1}(f.object)

testargs(f::FieldGradient,x::Point) = testargs(f.object,x)
return_value(f::FieldGradient,x::Point) = evaluate(f,testargs(f,x)...)
return_cache(f::FieldGradient,x::Point) = nothing
evaluate!(cache,f::FieldGradient,x::Point) = @abstractmethod
testvalue(::Type{FieldGradient{N,T}}) where {N,T} = FieldGradient{N}(testvalue(T))

# Default methods for arrays of points

function testargs(f::Field,x::AbstractArray{<:Point})
  y = copy(x)
  broadcast!(xi->first(testargs(f,xi)),y,x)
  (y,)
end

function return_cache(f::Field,x::AbstractArray{<:Point})
  T = return_type(f,testitem(x))
  s = size(x)
  ab = zeros(T,s)
  cb = CachedArray(ab)
  cf = return_cache(f,testitem(x))
  cb, cf
end

function evaluate!(c,f::Field,x::AbstractArray{<:Point})
  cb, cf = c
  sx = size(x)
  setsize!(cb,sx)
  r = cb.array
  for i in eachindex(x)
    @inbounds r[i] = evaluate!(cf,f,x[i])
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

#Field(f) = GenericField(f)
GenericField(f::Field) = f

testargs(a::GenericField,x::Point) = testargs(a.object,x)
return_value(a::GenericField,x::Point) = return_value(a.object,x)
return_cache(a::GenericField,x::Point) = return_cache(a.object,x)
evaluate!(cache,a::GenericField,x::Point) = evaluate!(cache,a.object,x)

function return_cache(f::FieldGradient{N,<:GenericField},x::Point) where N
  return_cache(FieldGradient{N}(f.object.object),x)
end

function evaluate!(c,f::FieldGradient{N,<:GenericField},x::Point) where N
  evaluate!(c,FieldGradient{N}(f.object.object),x)
end

# Make Field behave like a collection

Base.length(::Field) = 1
Base.size(::Field) = ()
Base.axes(::Field) = ()
Base.IteratorSize(::Type{<:Field}) = Base.HasShape{0}()
Base.eltype(::Type{T}) where T<:Field = T
Base.iterate(a::Field) = (a,nothing)
Base.iterate(a::Field,::Nothing) = nothing
Base.getindex(a::Field,i::Integer) =  (@check i == 1; a)
testitem(a::Field) = a

# Zero field

Base.zero(a::Field) = ZeroField(a)

"""
It represents `0.0*f` for a field `f`.
"""
struct ZeroField{F} <: Field
  field::F
end

return_cache(z::ZeroField,x::Point) = zero(return_type(z.field,x))
evaluate!(cache,z::ZeroField,x::Point) = cache
testvalue(::Type{ZeroField{F}}) where F = ZeroField(testvalue(F))

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

gradient(z::ZeroField) = ZeroField(gradient(z.field))

# Make Number behave like Field
#

# Number itself does not implement the Field interface since Number objects are not callable in Julia.
# Thus, wrapping a Number in a GenericField does not makes sense since the wrapped object
# is assumed to implement the Field interface.
# I think it is conceptually better to have ConstantField as struct otherwise we break the invariant
# "for any object wrapped in a GenericField we can assume that it implements the Field interface"
struct ConstantField{T<:Number} <: Field
  value::T
end

constant_field(a) = ConstantField(a)

Base.zero(::Type{ConstantField{T}}) where T = ConstantField(zero(T))

function evaluate!(c,f::ConstantField,x::Point)
  f.value
end

function return_cache(f::ConstantField,x::AbstractArray{<:Point})
  nx = size(x)
  c = fill(f.value,nx)
  CachedArray(c)
end

function evaluate!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = size(x)
  # This optimization is a bug if we include several ConstantField with different states
  # in the same array and we try to reuse cache between them.
  #if size(c) != nx
    setsize!(c,nx)
    fill!(c.array,f.value)
  #end
  c.array
end

function return_cache(f::FieldGradient{N,<:ConstantField},x::Point) where N
  gradient(f.object.value,Val(N))(x)
end

evaluate!(c,f::FieldGradient{N,<:ConstantField},x::Point) where N = c

function return_cache(f::FieldGradient{N,<:ConstantField},x::AbstractArray{<:Point}) where N
  CachedArray(gradient(f.object.value,Val(N)).(x))
end

function evaluate!(c,f::FieldGradient{N,<:ConstantField},x::AbstractArray{<:Point}) where N
  nx = size(x)
  if size(c) != nx
    setsize!(c,nx)
    fill!(c.array,zero(eltype(c)))
  end
  c.array
end

function lazy_map(::Operation{typeof(inv)},a::LazyArray{<:Fill{typeof(constant_field)}})
  v = a.args[1]
  vinv = lazy_map(inv,v)
  lazy_map(constant_field,vinv)
end

## Make Function behave like Field

return_cache(f::FieldGradient{N,<:Function},x::Point) where N = gradient(f.object,Val(N))
evaluate!(c,f::FieldGradient{N,<:Function},x::Point) where N = c(x)

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

function evaluate!(cache,c::OperationField,x::Point)
  ck, cf = cache
  lx = map((ci,fi) -> evaluate!(ci,fi,x),cf,c.fields)
  evaluate!(ck,c.op,lx...)
end

function return_value(c::OperationField,x::AbstractArray{<:Point})
  fx = map(f -> return_value(f,x),c.fields)
  c.op.(fx...)
end

function return_cache(c::OperationField,x::AbstractArray{<:Point})
  cf = map(fi -> return_cache(fi,x),c.fields)
  lx = map((ci,fi) -> evaluate!(ci,fi,x),cf,c.fields)
  ck = return_cache(c.op,map(testitem,lx)...)
  r = c.op.(lx...)
  ca = CachedArray(r)
  ca, ck, cf
end

function evaluate!(cache,c::OperationField,x::AbstractArray{<:Point})
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

evaluate!(cache,op::Operation,x::Field...) = OperationField(op.op,x)
return_value(op::Broadcasting{<:Operation},x::Field...) = OperationField(op.f.op,x)
evaluate!(cache,op::Broadcasting{<:Operation},x::Field...) = OperationField(op.f.op,x)

# Define some well known operations

for op in (:+,:-,:*,:/,:⋅,:⊙,:⊗,:inv,:det,:meas,:pinvJt,:tr,:grad2curl,:symmetric_part,:transpose)
  @eval ($op)(a::Field...) = Operation($op)(a...)
end

transpose(f::Field) = f

for op in (:+,:-,:*,:/,:⋅,:⊙,:⊗)
  @eval ($op)(a::Field,b::Number) = Operation($op)(a,ConstantField(b))
  @eval ($op)(a::Number,b::Field) = Operation($op)(ConstantField(a),b)
end

#*(A::Number, B::Field) = ConstantField(A)*B
#*(A::Field, B::Number) = A*ConstantField(B)
#⋅(A::Number, B::Field) = ConstantField(A)⋅B
#⋅(A::Field, B::Number) = A⋅ConstantField(B)

#*(A::Function, B::Field) = GenericField(A)*B
#*(A::Field, B::Function) = GenericField(B)*A

# Gradient of the sum
for op in (:+,:-)
  @eval begin
    function gradient(a::OperationField{typeof($op)})
      f = a.fields
      g = map( gradient, f)
      $op(g...)
    end
  end
end

# Gradient of the product

function product_rule(fun,f1,f2,∇f1,∇f2)
  msg = "Product rule not implemented for product $fun between types $(typeof(f1)) and $(typeof(f2))"
  @notimplemented msg
end

function product_rule(::typeof(*),f1::Real,f2::Real,∇f1,∇f2)
  ∇f1*f2 + f1*∇f2
end

for op in (:*,:⋅)
  @eval begin
     function product_rule(::typeof($op),f1::Real,f2::VectorValue,∇f1,∇f2)
       ∇f1⊗f2 + ∇f2*f1
     end

     function product_rule(::typeof($op),f1::VectorValue,f2::Real,∇f1,∇f2)
       product_rule(*,f2,f1,∇f2,∇f1)
     end
  end
end

function product_rule(::typeof(⋅),f1::VectorValue,f2::VectorValue,∇f1,∇f2)
  ∇f1⋅f2 + ∇f2⋅f1
end

function product_rule(::typeof(⋅),f1::TensorValue,f2::VectorValue,∇f1,∇f2)
  ∇f1⋅f2 + ∇f2⋅transpose(f1)
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
Base.:∘(f::Field,g::Field) = Operation(f)(g)
evaluate!(cache,::Broadcasting{typeof(∘)},f::Field,g::Field) = f∘g

# Integration

"""
Integration of a given field in the "physical" space
"""
function integrate(a::Field,x::AbstractVector{<:Point},w::AbstractVector{<:Real})
  cache = return_cache(integrate,a,x,w)
  evaluate!(cache,integrate,a,x,w)
end

"""
Integration of a given field in the "reference" space
"""
function integrate(a::Field,q::AbstractVector{<:Point},w::AbstractVector{<:Real},j::Field)
  cache = return_cache(integrate,a,q,w,j)
  evaluate!(cache,integrate,a,q,w,j)
end

function return_cache(::typeof(integrate),a,x,w)
  ca = return_cache(a,x)
  ax = return_value(a,x)
  ck = return_cache(IntegrationMap(),ax,w)
  ca, ck
end

function evaluate!(cache,::typeof(integrate),a,x,w)
  ca, ck = cache
  ax = evaluate!(ca,a,x)
  evaluate!(ck,IntegrationMap(),ax,w)
end

function return_cache(::typeof(integrate),a,q,w,j)
  ca = return_cache(a,q)
  cj = return_cache(j,q)
  aq = return_value(a,q)
  jq = return_value(j,q)
  ck = return_cache(IntegrationMap(),aq,w,jq)
  ca, cj, ck
end

function evaluate!(cache,::typeof(integrate),a,q,w,j)
  ca, cj, ck = cache
  aq = evaluate!(ca,a,q)
  jq = evaluate!(cj,j,q)
  evaluate!(ck,IntegrationMap(),aq,w,jq)
end

struct IntegrationMap <: Map end

function evaluate!(cache,k::IntegrationMap,ax::AbstractVector,w)
  T = typeof( testitem(ax)*testitem(w) + testitem(ax)*testitem(w) )
  z = zero(T)
  r = z
  @check length(ax) == length(w)
  @inbounds for i in eachindex(ax)
    r += ax[i]*w[i]
  end
  r
end

function evaluate!(cache,k::IntegrationMap,aq::AbstractVector,w,jq::AbstractVector)
  T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq)) + testitem(aq)*testitem(w)*meas(testitem(jq)) )
  z = zero(T)
  @check length(aq) == length(w)
  @check length(aq) == length(jq)
  @inbounds for i in eachindex(aq)
    z += aq[i]*w[i]*meas(jq[i])
  end
  z
end

function return_cache(k::IntegrationMap,ax::AbstractArray,w)
  T = typeof( testitem(ax)*testitem(w) + testitem(ax)*testitem(w) )
  r = zeros(T,size(ax)[2:end])
  CachedArray(r)
end

function evaluate!(cache,k::IntegrationMap,ax::AbstractArray,w)
  setsize!(cache,size(ax)[2:end])
  r = cache.array
  @check size(ax,1) == length(w)
  @inbounds for j in CartesianIndices(r)
    rj = zero(eltype(r))
    for p in 1:length(w)
      rj += ax[p,j]*w[p]
    end
    r[j] = rj
  end
  r
end

function return_value(k::IntegrationMap,aq::AbstractArray,w,jq::AbstractVector)
  if size(aq,1) == length(w) && size(aq,1) == length(jq)
    evaluate(k,aq,w,jq)
  else
    c = return_cache(k,aq,w,jq)
    c.array
  end
end

function return_cache(k::IntegrationMap,aq::AbstractArray,w,jq::AbstractVector)
  T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq)) + testitem(aq)*testitem(w)*meas(testitem(jq)) )
  r = zeros(T,size(aq)[2:end])
  CachedArray(r)
end

function evaluate!(cache,k::IntegrationMap,aq::AbstractArray,w,jq::AbstractVector)
  setsize!(cache,size(aq)[2:end])
  r = cache.array
  @check size(aq,1) == length(w) || size(aq,1) == 0
  @check size(aq,1) == length(jq) || size(aq,1) == 0
  fill!(r,zero(eltype(r)))
  cis = CartesianIndices(r)
  @inbounds for p in 1:length(w)
    dV = meas(jq[p])*w[p]
    for j in cis
      r[j] += aq[p,j]*dV
    end
  end
  r
end

function return_value(k::IntegrationMap,aq::AbstractArray{S,3} where S,w,jq::AbstractVector)
  T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq)) + testitem(aq)*testitem(w)*meas(testitem(jq)) )
  r = zeros(T,size(aq)[2:end])
  r
end

function return_cache(k::IntegrationMap,aq::AbstractArray{S,3} where S,w,jq::AbstractVector)
  T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq)) + testitem(aq)*testitem(w)*meas(testitem(jq)) )
  r = zeros(T,size(aq)[2:end])
  s = zeros(typeof(meas(testitem(jq))),length(jq))
  CachedArray(r), CachedArray(s)
end

function evaluate!(cache,k::IntegrationMap,aq::AbstractArray{S,3} where S, w,jq::AbstractVector)
  cache_r, cache_s = cache
  np, ni, nj = size(aq)
  setsize!(cache_r,(ni,nj))
  setsize!(cache_s,(np,))
  r = cache_r.array
  dV = cache_s.array
  @check np == length(w) || np == 0
  @check np == length(jq) || np == 0
  @inbounds for p in 1:np
    dV[p] = meas(jq[p])*w[p]
  end
  #fill!(r,zero(eltype(r)))
  @inbounds for j in 1:nj
    for i in 1:ni
      rij = zero(eltype(aq))
      for p in 1:np
        rij += aq[p,i,j]*dV[p]
      end
      r[i,j] = rij
    end
  end
  r
end

function evaluate!(cache,k::IntegrationMap,aq::AbstractMatrix, w,jq::AbstractVector)
  np, ni = size(aq)
  setsize!(cache,(ni,))
  r = cache.array
  @check np == length(w) || np == 0
  @check np == length(jq) || np == 0
  fill!(r,zero(eltype(r)))
  @inbounds for p in 1:np
    dV = meas(jq[p])*w[p]
    for i in 1:ni
      r[i] += aq[p,i]*dV
    end
  end
  r
end

# Testers

"""
    test_field(
      f::Union{Field,AbstractArray{<:Field}},
      x,
      v,
      cmp=(==);
      grad=nothing,
      gradgrad=nothing)

Function used to test the field interface. `v` is an array containing the expected
result of evaluating the field `f` at the point or vector of points `x`. The comparison is performed using
the `cmp` function. For fields objects that support the `gradient` function, the keyword
argument `grad` can be used. It should contain the result of evaluating `gradient(f)` at x.
Idem for `gradgrad`. The checks are performed with the `@test` macro.
"""
function test_field(f::Field, x, v, cmp=(==); grad=nothing, gradgrad=nothing)
  test_map(v,f,x;cmp=cmp)
  if grad != nothing
    test_map(grad,∇(f),x;cmp=cmp)
  end
  if gradgrad != nothing
    test_map(gradgrad,∇∇(f),x;cmp=cmp)
  end
  true
end

struct VoidFieldMap <: Map
  isvoid::Bool
end

Arrays.evaluate!(cache,k::VoidFieldMap,b) = VoidField(b,k.isvoid)

struct VoidField{F} <: Field
  field::F
  isvoid::Bool
end

function return_cache(f::VoidField,x::Point)
  c = return_cache(f.field,x)
  fx = evaluate!(c,f.field,x)
  c, zero(fx)
end

function evaluate!(cache,f::VoidField,x::Point)
  c,z = cache
  if f.isvoid
    z
  else
    evaluate!(c,f.field,x)
  end
end

function return_cache(f::VoidField,x::AbstractVector{<:Point})
  c = return_cache(f.field,x)
  fx = evaluate!(c,f.field,x)
  z = similar(fx)
  c, CachedArray(z)
end

function evaluate!(cache,f::VoidField,x::AbstractVector{<:Point})
  c,z = cache
  if f.isvoid
    setsize!(z,size(x))
    fill!(z,zero(eltype(z)))
    z.array
  else
    evaluate!(c,f.field,x)
  end
end

testvalue(::Type{VoidField{F}}) where F = VoidField(testvalue(F),false)

gradient(z::VoidField) = VoidField(gradient(z.field),z.isvoid)

function lazy_map(::typeof(evaluate),a::LazyArray{<:Fill{VoidFieldMap}},x::AbstractArray)
  p = a.maps.value
  @notimplementedif p.isvoid
  lazy_map(evaluate,a.args[1],x)
end

struct VoidBasisMap <: Map
  isvoid::Bool
end

Arrays.evaluate!(cache,k::VoidBasisMap,b) = VoidBasis(b,k.isvoid)

struct VoidBasis{T,N,A} <: AbstractArray{T,N}
  basis::A
  isvoid::Bool
  function VoidBasis(basis::AbstractArray{T,N},isvoid::Bool) where {T,N}
    new{T,N,typeof(basis)}(basis,isvoid)
  end
end

function Base.size(a::VoidBasis)
  if a.isvoid
    0 .* size(a.basis)
  else
    size(a.basis)
  end
end

Base.IndexStyle(::Type{<:VoidBasis}) = IndexLinear()

function Base.getindex(a::VoidBasis,i::Integer)
  if a.isvoid
    @unreachable "Unable to access 0-length array"
  else
    a.basis[i]
  end
end

Arrays.testitem(a::VoidBasis) = testitem(a.basis)

function _zero_size(a::VoidBasis{T,1} where T)
  (0,)
end

function _zero_size(a::VoidBasis{T,2} where T)
  @check size(a,1) in (0,1)
  (1,0)
end

function Fields.return_cache(a::VoidBasis,x::Point)
  cb = return_cache(a.basis,x)
  bx = return_value(a.basis,x)
  zs = _zero_size(a)
  r = similar(bx,zs)
  cb,r
end

function Fields.return_cache(a::VoidBasis,x::Field)
  @notimplementedif ndims(a) != 1
  cb = return_cache(a.basis,x)
  bx = return_value(a.basis,x)
  r = similar(bx,(0,))
  cb,r
end

function Fields.return_cache(a::VoidBasis,x::AbstractVector{<:Point})
  cb = return_cache(a.basis,x)
  bx = return_value(a.basis,x)
  zs = _zero_size(a)
  r = similar(bx,(length(x),zs...))
  cb,r
end

function Fields.return_cache(a::VoidBasis,v::AbstractVector{<:Field})
  @notimplementedif ndims(a) != 1
  cb = return_cache(a.basis,v)
  bx = return_value(a.basis,v)
  r = similar(bx,(0,length(v)))
  cb,r
end

for T in (:Point,:Field,:(AbstractVector{<:Point}),:(AbstractVector{<:Field}))
  @eval begin

    function Fields.evaluate!(cache,a::VoidBasis,x::$T)
      cb, r = cache
      if a.isvoid
        r
      else
        evaluate!(cb,a.basis,x)
      end
    end

  end
end

function Fields.evaluate!(cache,k::Broadcasting{typeof(∇)},a::VoidBasis)
  VoidBasis(k(a.basis),a.isvoid)
end

function Fields.evaluate!(cache,k::Broadcasting{typeof(∇∇)},a::VoidBasis)
  VoidBasis(k(a.basis),a.isvoid)
end

function lazy_map(::typeof(evaluate),a::LazyArray{<:Fill{VoidBasisMap}},x::AbstractArray)
  p = a.maps.value
  @notimplementedif p.isvoid
  lazy_map(evaluate,a.args[1],x)
end
