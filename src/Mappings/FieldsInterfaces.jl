const Point{D,T} = VectorValue{D,T}

abstract type Field <: Mapping end

const FieldArray{T,N} = AbstractArray{T,N} where {T<:Field,N}

const FieldOrFieldArray = Union{Field,FieldArray}

function return_cache(f::Field,x::AbstractArray{<:Point})
  T = return_type(f,first(x))
  s = size(x)
  ab = zeros(T,s)
  cb = CachedArray(ab)
  cf = return_cache(f,first(x))
  cb, cf
end

function evaluate!(c,f::Field,x::AbstractArray{<:Point})
  cb, cf = c
  sx = size(x)
  setsize!(cb,sx)
  for i in eachindex(x)
    @inbounds cb[i] = evaluate!(cf,f,x[i])
  end
  cb.array
end

# Differentiation

function gradient end

const ∇ = gradient

@inline function evaluate_gradient!(cache,f,x)
  @abstractmethod
end

@inline function return_gradient_cache(f,x)
  nothing
end

@inline function evaluate_hessian!(cache,f,x)
  @abstractmethod
end

@inline function return_hessian_cache(f,x)
  nothing
end

@inline function return_gradient_type(::Type{T},x::Point) where T
  typeof(outer(zero(x),zero(T)))
end

@inline function return_hessian_type(::Type{T},x::Point) where T
  typeof(outer(zero(x),zero(return_gradient_type(T,x))))
end

# GenericField

struct GenericField{T} <: Field
  object::T
end

@inline Field(f) = GenericField(f)
@inline GenericField(f::Field) = f

@inline return_type(::Type{<:GenericField},::Type{T}) where T<:Field = T
@inline return_type(::Type{<:GenericField},::Type{T}) where T = GenericField{T}

@inline return_cache(a::GenericField,x::Point) = return_cache(a.object,x)
@inline return_cache(a::GenericField,x::AbstractArray{<:Point}) = return_cache(a.object,x)

@inline return_type(a::GenericField,x::Point) = return_type(a.object,x)
@inline return_type(a::GenericField,x::AbstractArray{<:Point}) = return_type(a.object,x)

@inline evaluate!(cache,a::GenericField,x::Point) = evaluate!(cache,a.object,x)
@inline evaluate!(cache,a::GenericField,x::AbstractArray{<:Point}) = evaluate!(cache,a.object,x)

# Make Field behave like a collection

@inline Base.length(::Field) = 1
@inline Base.size(::Field) = ()
@inline Base.axes(::Field) = ()
@inline Base.IteratorSize(::Type{<:Field}) = Base.HasShape{0}()
@inline Base.eltype(::Type{T}) where T<:Field = T
@inline Base.iterate(a::Field) = (a,nothing)
@inline Base.iterate(a::Field,::Nothing) = nothing

# Zero field

B@inline ase.zero(a::Field) = ZeroField(a)

struct ZeroField{F} <: Field
  field::F
end

@inline return_type(z::ZeroField,x::Point) = return_type(z.field,x)
@inline return_type(z::ZeroField,x::AbstractArray{<:Point}) = return_type(z.field,x)

@inline return_cache(z::ZeroField,x::Point) = zero(return_type(z.field,x))

function return_cache(z::ZeroField,x::AbstractArray{<:Point})
  E = return_type(z.field,first(x))
  c = zeros(E,length(x))
  CachedArray(c)
end

@inline evaluate!(cache,z::ZeroField,x::Point) = cache

function evaluate!(c,f::ZeroField,x::AbstractArray{<:Point})
  nx = length(x)
  if size(c) != nx
    setsize!(c,(nx,))
    c .= zero(eltype(c))
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

const ConstantField{T} = GenericField{T} where T<:Number

@inline function evaluate!(c,f::ConstantField,x::Point)
  f.object
end

function return_type(f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  c = zeros(typeof(f.object),nx)
  typeof(c)
end

function return_cache(f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  c = zeros(typeof(f.object),nx)
  CachedArray(c)
end

function evaluate!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  setsize!(c,(nx,))
  r = c.array
  for i in eachindex(x)
    @inbounds r[i] = f.object
  end
  r
end

@inline function return_gradient_cache(f::ConstantField,x::Point)
  gradient(f.object)(x)
end

@inline function return_gradient_cache(f::ConstantField,x::AbstractArray{<:Point})
  CachedArray(gradient(f.object).(x))
end

@inline evaluate_gradient!(c,f::ConstantField,x::Point) = c

function evaluate_gradient!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  if size(c) != nx
    setsize!(c,(nx,))
    c .= zero(eltype(c))
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
  nx = length(x)
  if size(c) != nx
    setsize!(c,(nx,))
    c .= zero(eltype(c))
  end
  c.array
end

# Make Function behave like Field

const FunctionField{F} = GenericField{F} where F<:Function

function return_cache(f::FunctionField,x::AbstractArray{<:Point})
  nx = length(x)
  Te = eltype(x)
  c = zeros(return_type(f.object,Te),nx)
  CachedArray(c)
end

function evaluate!(c,f::FunctionField,x::AbstractArray{<:Point})
  nx = length(x)
  setsize!(c,(nx,))
  for i in eachindex(x)
    c[i] = f.object(x[i])
  end
  c.array
end

@inline function return_gradient_cache(f::FunctionField,x::Point)
  gradient(f.object)

end

function return_gradient_cache(f::FunctionField,x::AbstractArray{<:Point})
  gf = gradient(f.object)
  nx = length(x)
  Te = eltype(x)
  c = zeros(return_type(gf,Te),nx)
  # gf, CachedArray(c)
  gf, CachedArray(c)
end

@inline evaluate_gradient!(c,f::FunctionField,x::Point) = c(x)

function evaluate_gradient!(cache,f::FunctionField,x::AbstractArray{<:Point})
  gf, c = cache
  nx = length(x)
  setsize!(c,(nx,))
  for i in eachindex(x)
    c[i] = gf(x[i])
  end
  c.array
end

# Differentiation

struct FieldGradient{F} <: Field
  object::F
end

@inline gradient(f::Field) = FieldGradient(f)

@inline gradient(f::GenericField{FieldGradient}) = FieldHessian(f.object.object)

@inline evaluate!(cache,f::FieldGradient,x::Point) = evaluate_gradient!(cache,f.object,x)
@inline evaluate!(cache,f::FieldGradient,x::AbstractArray{<:Point}) = evaluate_gradient!(cache,f.object,x)

@inline return_cache(f::FieldGradient,x::Point) = return_gradient_cache(f.object,x)
@inline return_cache(f::FieldGradient,x::AbstractArray{<:Point}) = return_gradient_cache(f.object,x)

struct FieldHessian{F} <: Field
  object::F
end

gradient(f::FieldGradient) = FieldHessian(f.object)

@inline evaluate!(cache,f::FieldHessian,x::Point) = evaluate_hessian!(cache,f.object,x)
@inline evaluate!(cache,f::FieldHessian,x::AbstractArray{<:Point}) = evaluate_hessian!(cache,f.object,x)

@inline return_cache(f::FieldHessian,x::Point) = return_hessian_cache(f.object,x)
@inline return_cache(f::FieldHessian,x::AbstractArray{<:Point}) = return_hessian_cache(f.object,x)

@inline function gradient(f::FieldHessian)
  @unreachable "Default implementation of 3rt order derivatives not available"
end

# Operations

struct OperationField{O,F} <: Field
  op::O
  fields::F
end

function return_type(c::OperationField,x::Point)
  return_type(c.op,evaluate(c.fields,x)...)
  # return_type(c.op,map(typeof,evaluate(c.fields,x))...)
end

function return_cache(c::OperationField,x::Point)
  cl = return_caches(c.fields,x)
  lx = evaluate!(cl,c.fields,x)
  ck = return_cache(c.op,lx)
  ck, cl
end

function return_cache(c::OperationField,x::AbstractArray{<:Point})
  cl = return_caches(c.fields,x)
  lx = evaluate!(cl,c.fields,x)
  ck = CachedArray(zero(c.op.(lx...)))
  ck, cl
end

@inline function evaluate!(cache,c::OperationField,x::AbstractArray{<:Point})
  ck, cf = cache
  sx = size(x)
  setsize!(ck,sx)
  lx = evaluate!(cf,c.fields,x)
  for i in eachindex(x)
    @inbounds ck.array[i] = c.op(map(lxi -> lxi[i], lx)...)
  end
  ck.array
end

@inline function evaluate!(cache,c::OperationField,x::Point)
  ck, cf = cache
  lx = evaluate!(cf,c.fields,x)
  # evaluate!(ck,c.op,lx)
  c.op(lx...)
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
@inline *(A::Number, B::Field) = GenericField(A)*B
@inline *(A::Field, B::Number) = GenericField(B)*A
@inline *(A::Function, B::Field) = GenericField(A)*B
@inline *(A::Field, B::Function) = GenericField(B)*A

function gradient(a::OperationField{typeof(⋅)})
  f = a.fields
  if length(f) != 2 @notimplemented end
  f1, f2 = f
  g1, g2 = map(gradient, f)
  g1⋅f2+f1⋅g2
end

function gradient(a::OperationField{typeof(*)})
  f = a.fields
  if length(f) != 2 @notimplemented end
  f1, f2 = f
  g1, g2 = map(gradient, f)
  g1⋅f2+f1⋅g2
end

# Chain rule
function gradient(f::OperationField{<:Field})
  a = f.op
  @notimplementedif length(f.fields) != 1
  b, = f.fields
  _x = ∇(a)∘b
  _y = ∇(b)
  _x⋅_y
end

# Composition

@inline Base.:∘(f::Field,g::Field) = Operation(f)(g)

# Other operations

@inline transpose(f::Field) = f
@inline Base.copy(f::Field) = f
@inline *(f::Field,g::Field) = Operation(*)(f,g)#⋅g

# Testers

function test_field(
  f::FieldOrFieldArray,
  x::Tuple,
  v,cmp=(==);
  grad=nothing,
  hessian=nothing)

  x, = x

  @test isa(x,Union{Point,AbstractArray{<:Point}})

  w = evaluate(f,x)

  @test cmp(w,v)
  @test typeof(w) == return_type(f,x)

  cf = return_cache(f,x)
  r = evaluate!(cf,f,x)
  @test cmp(r,v)

  if x isa AbstractArray{<:Point}

    _x = vcat(x,x)
    _v = vcat(v,v)
    _w = evaluate!(cf,f,_x)
    @test cmp(_w,_v)
  end

  if isa(f,Field)
    test_mapping(f,(x,),v,cmp)
  end

  if grad != nothing
    g = gradient(f)
    if typeof(f) <: Field
      @test g isa Field
    elseif typeof(f) <: AbstractArray{<:Field}
      @test g isa AbstractArray{<:Field}
    end
    test_field(g,(x,),grad,cmp,grad=hessian)
  end
end

@inline function test_field(f::FieldOrFieldArray,x,v,cmp=(==);grad=nothing,hessian=nothing)
  test_field(f,(x,),v,cmp;grad=grad,hessian=hessian)
end
