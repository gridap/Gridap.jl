# Meeting 22 Sep
# integrate(f::Field,x::AbstractVector{<:Point},w::AbstractVector{<:Real}) = sum( f(x) .* w  )
# integrate(f::AbstractArray{<:Field},x::AbstractVector{<:Point},w::AbstractVector{<:Real}) =

const Point{D,T} = VectorValue{D,T}

abstract type Field <: Mapping end

const FieldArray{T,N} = AbstractArray{T,N} where {T<:Field,N}

const FieldOrFieldArray = Union{Field,FieldArray}

function evaluate!(cache,f::Field,x)
  @abstractmethod
end

# Differentiation

function gradient end

const ∇ = gradient

function evaluate_gradient!(cache,f,x)
  @abstractmethod
end

function return_gradient_cache(f,x)
  nothing
end

function evaluate_hessian!(cache,f,x)
  @abstractmethod
end

function return_hessian_cache(f,x)
  nothing
end

function return_gradient_type(::Type{T},x::Point) where T
  typeof(outer(zero(x),zero(T)))
end

function return_hessian_type(::Type{T},x::Point) where T
  typeof(outer(zero(x),zero(return_gradient_type(T,x))))
end

# function return_type(f::Field,x::Point)
#   typeof(evaluate(f,x))
# end

# GenericField

struct GenericField{T} <: Field
  object::T
end

Field(f) = GenericField(f)
GenericField(f::Field) = f

@inline return_type(::Type{<:GenericField},::Type{T}) where T<:Field = T
@inline return_type(::Type{<:GenericField},::Type{T}) where T = GenericField{T}

@inline return_cache(a::GenericField,x::Point) = return_cache(a.object,x)
@inline return_cache(a::GenericField,x::AbstractArray{<:Point}) = return_cache(a.object,x)

@inline return_type(a::GenericField,x::Point) = return_type(a.object,x)
@inline return_type(a::GenericField,x::AbstractArray{<:Point}) = return_type(a.object,x)

@inline evaluate!(cache,a::GenericField,x::Point) = evaluate!(cache,a.object,x)
@inline evaluate!(cache,a::GenericField,x::AbstractArray{<:Point}) = evaluate!(cache,a.object,x)

# Make Field behave like a collection

Base.length(::Field) = 1
Base.size(::Field) = ()
Base.axes(::Field) = ()
Base.IteratorSize(::Type{<:Field}) = Base.HasShape{0}()
Base.eltype(::Type{T}) where T<:Field = T
Base.iterate(a::Field) = (a,nothing)
Base.iterate(a::Field,::Nothing) = nothing

# Zero field

Base.zero(a::Field) = ZeroField(a)

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

gradient(z::ZeroField) = ZeroField(gradient(z.field))

# function gradient(f::BroadcastField{GenericField{<:OperationMapping{<:Field}}})#,Tuple{<:Field}}})
#   a = f.object.k
#   @assert length(f.object.l) == 1
#   b, = f.object.l
#   _x = ∇(a)∘b
#   _y = ∇(b)
#   # @santiagobadia : How can we express this product ????
# end

# function gradient(a::GenericField{<:OperationMapping{typeof(*)}})
# function gradient(f::BroadcastField{GenericField})#{OperationMapping{BroadcastFunction{typeof(*)}}})
#   f = a.object.l
#   if length(f) != 2 @notimplemented end
#   f1, f2 = f
#   g1, g2 = map( gradient, f)
#   g1.⋅f2+f1.⋅g2
# end

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

function return_gradient_cache(f::ConstantField,x::Point)
  gradient(f.object)(x)
end

function return_gradient_cache(f::ConstantField,x::AbstractArray{<:Point})
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

function return_hessian_cache(f::ConstantField,x::Point)
  hessian(f.object)(x)
end

function return_hessian_cache(f::ConstantField,x::AbstractArray{<:Point})
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

function_field(f::Function) =  GenericField(f)

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

# function evaluate!(cache,a::Function,x::AbstractArray{<:Point})
#   evaluate!(cache,FunctionField(a),x)
# end

# function return_cache(f::Function,x::AbstractArray{<:Point})
#   return_cache(FunctionField(f),x)
# end

function return_gradient_cache(f::FunctionField,x::Point)
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

struct Gradient{F} <: Field
  object::F
end

gradient(f::Field) = Gradient(f)

gradient(f::GenericField{Gradient}) = Hessian(f.object.object)

@inline evaluate!(cache,f::Gradient,x::Point) = evaluate_gradient!(cache,f.object,x)
@inline evaluate!(cache,f::Gradient,x::AbstractArray{<:Point}) = evaluate_gradient!(cache,f.object,x)

@inline return_cache(f::Gradient,x::Point) = return_gradient_cache(f.object,x)
@inline return_cache(f::Gradient,x::AbstractArray{<:Point}) = return_gradient_cache(f.object,x)

struct Hessian{F} <: Field
  object::F
end

gradient(f::Gradient) = Hessian(f.object)

@inline evaluate!(cache,f::Hessian,x::Point) = evaluate_hessian!(cache,f.object,x)
@inline evaluate!(cache,f::Hessian,x::AbstractArray{<:Point}) = evaluate_hessian!(cache,f.object,x)

@inline return_cache(f::Hessian,x::Point) = return_hessian_cache(f.object,x)
@inline return_cache(f::Hessian,x::AbstractArray{<:Point}) = return_hessian_cache(f.object,x)

@inline function gradient(f::Hessian)
  @unreachable "Default implementation of 3rt order derivatives not available"
end

# Operations

struct OperationField{O,F} <: Field
  op::O
  fields::F
end

function return_type(c::OperationField,x::Point)
  return_type(c.op,map(typeof,evaluate(c.fields,x))...)
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
    @inbounds ck.array[i] = c.op(map(lxi -> lxi[i], lx)...) #evaluate!(cfi,f[i],x)
  end
  # Misssing result cache
  # c.op.(lx...)
  ck.array
end

@inline function evaluate!(cache,c::OperationField,x::Point)
  ck, cf = cache
  lx = evaluate!(cf,c.fields,x)
  # evaluate!(ck,c.op,lx)
  c.op(lx...)
end

evaluate!(cache,op::Operation,x::Field...) = OperationField(op.op,x)

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
*(A::Number, B::Field) = GenericField(A)*B
*(A::Field, B::Number) = GenericField(B)*A
*(A::Function, B::Field) = GenericField(A)*B
*(A::Field, B::Function) = GenericField(B)*A

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

Base.:∘(f::Field,g::Field) = Operation(f)(g)

# Broadcast fields

# struct BroadcastField{T<:Field} <: Field
#   field::T
# end

# BroadcastField(f::BroadcastField) = f

function return_cache(f::Field,x::AbstractArray{<:Point})
  T = return_type(f,first(x))
  s = size(x)
  ab = zeros(T,s)
  cb = CachedArray(ab)
  cf = return_cache(f,first(x))
  cb, cf
end

# return_type(f::BroadcastField,x::Point) = return_type(f.field,x)

# function return_cache(f::BroadcastField,x::Point)
#   return_cache(f.field,x)
# end

function evaluate!(c,f::Field,x::AbstractArray{<:Point})
  cb, cf = c
  sx = size(x)
  setsize!(cb,sx)
  for i in eachindex(x)
    @inbounds cb[i] = evaluate!(cf,f,x[i])
  end
  cb.array
end

# evaluate!(c,f::BroadcastField,x::Point) = evaluate!(c,f.field,x)

# gradient(b::BroadcastField) = BroadcastField(gradient(b.field))

# for op in (:+,:-,:*,:⋅,:inv,:det)
#   # @eval ($op)(a::BroadcastField...) = BroadcastField(Operation($op)(map( f -> f.field,a)...))
#   @eval ($op)(a::BroadcastField...) = BroadcastField(Operation(BroadcastFunction($op))(map( f -> f.field,a)...))
# end

# struct BroadcastFunction{T}
#   op::T
# end

# (f::BroadcastFunction)(x) = f.op.(x)

# evaluate!(cache,f::BroadcastFunction,x...) = f.op.(x...)

# Other operations

transpose(f::Field) = f
Base.copy(f::Field) = f
(*)(f::Field,g::Field) = f⋅g

# @santiagobadia : Just to make things work for the moment
(*)(f::VectorValue,g::VectorValue) = f⋅g
(*)(f::TensorValue,g::TensorValue) = f⊙g

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
    #  np, = size(w)
    #  @test size(x) == size(w)

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

function test_field(f::FieldOrFieldArray,x,v,cmp=(==);grad=nothing,hessian=nothing)
  test_field(f,(x,),v,cmp;grad=grad,hessian=hessian)
end
