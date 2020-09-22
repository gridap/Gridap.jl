# Meeting 22 Sep
# integrate(f::Field,x::AbstractVector{<:Point},w::AbstractVector{<:Real}) = sum( f(x) .* w  )
# integrate(f::AbstractArray{<:Field},x::AbstractVector{<:Point},w::AbstractVector{<:Real}) =
# linear_combination 17:21 per un tema d'optimitzacio 17:21 en el draft tinc transpose(i_to_f)*i_to_val
# linear_combination(i_to_f,i_to_val) = transpose(i_to_f)*i_to_val
# apply(linear_compination, cell_to_i_to_f,cell_to_i_to_val) es pot potmitizarA
# apply((a,b) ->transpose(a)*b, cell_to_i_to_f,cell_to_i_to_val) dona el mateix resultat, pero no es facil optimitzar-ho pq tens una anonymous function

const Point{D,T} = VectorValue{D,T}

abstract type Field <: Mapping end

const FieldVector{T} = AbstractVector{T} where T<:Field

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

# GenericField

struct GenericField{T} <: Field
  object::T
end

Field(f) = GenericField(f)
GenericField(f::Field) = f
return_type(::Type{<:GenericField},::Type{T}) where T<:Field = T
return_type(::Type{<:GenericField},::Type{T}) where T = GenericField{T}

@inline return_cache(a::GenericField,x) = return_cache(a.object,x)
@inline return_type(a::GenericField,x) = return_type(a.object,x)
@inline evaluate!(cache,a::GenericField,x) = evaluate!(cache,a.object,x)

# Make Field behave like a collection

# Base.length(::Field) = 1
# Base.size(::Field) = ()
# Base.axes(::Field) = ()
# Base.IteratorSize(::Type{<:Field}) = Base.HasShape{0}()
# Base.eltype(::Type{T}) where T<:Field = T
# Base.iterate(a::Field) = (a,nothing)
# Base.iterate(a::Field,::Nothing) = nothing

# Zero field

Base.zero(a::Field) = ZeroField(a)

struct ZeroField{F} <: Field
  field::F
end

@inline return_type(z::ZeroField,x) = return_type(z.field,x)
@inline evaluate!(cache,z::ZeroField,x::Point) = zero(return_type(z.field,x))

@inline function evaluate_gradient!(cache,z::ZeroField,x::Point)
 outer(zero(return_type(z.field)),zero(x))
end

gradient(z::ZeroField) = ZeroField(gradient(z.field))

# function evaluate!(cache,z::ZeroField,x::AbstractArray{<:Point})
#   T = return_type(z.f,testitem(x))
#   zeros(T,length(x)) # TODO cache
# end

# Operations

evaluate!(cache,op::Operation,x::Field...) = GenericField(OperationMapping(op.op,x))

for op in (:+,:-,:*,:⋅,:inv,:det)
  @eval ($op)(a::Field...) = Operation($op)(a...)
end

# Operation rules e.g.
for op in (:+,:-)
  @eval begin
    function gradient(a::GenericField{<:OperationMapping{typeof($op)}})
      f = a.object.l
      g = map( gradient, f)
      $op(g...)
    end
  end
end

function gradient(a::GenericField{<:OperationMapping{typeof(⋅)}})
  f = a.object.l
  if length(f) != 2 @notimplemented end
  f1, f2 = f
  g1, g2 = map( gradient, f)
  g1⋅f2+f1⋅g2
end

function gradient(a::GenericField{<:OperationMapping{typeof(*)}})
  f = a.object.l
  if length(f) != 2 @notimplemented end
  f1, f2 = f
  g1, g2 = map( gradient, f)
  g1⋅f2+f1⋅g2
end

# Chain rule
function gradient(f::GenericField{<:OperationMapping{<:Field,Tuple{<:Field}}})
  a = f.object.k
  @assert length(f.object.l) == 1
  b, = f.object.l
  (∇a∘b)⋅∇b
end

# Make Number behave like Field

const ConstantField{T} = GenericField{T} where T<:Number

@inline function evaluate!(c,f::ConstantField,x::Point)
  f.object
end

function return_gradient_cache(f::ConstantField,x)
  gradient(f.object)(x)
end

@inline evaluate_gradient!(c,f::ConstantField,x) = c

function return_hessian_cache(f::ConstantField,x)
  hessian(f.object)(x)
end

@inline evaluate_hessian!(c,f::ConstantField,x) = c

# function return_cache(f::ConstantField,x::AbstractArray{<:Point})
#   nx = length(x)
#   c = zeros(typeof(f.v),nx)
#   CachedArray(c)
# end

# function evaluate!(c,f::ConstantField,x::AbstractArray{<:Point})
#   nx = length(x)
#   setsize!(c,(nx,))
#   r = c.array
#   for i in eachindex(x)
#     @inbounds r[i] = f.v
#   end
#   r
# end

# function return_gradient_cache(f::ConstantField{T},x) where T
#   E = return_gradient_type(T,first(x))
#   c = zeros(E,length(x))
#   CachedArray(c)
# end

# function evaluate_gradient!(c,f::ConstantField,x)
#   nx = length(x)
#   if size(c) != nx
#     setsize!(c,(nx,))
#     c .= zero(eltype(c))
#   end
#   c
# end

# function return_hessian_cache(f::ConstantField{T},x) where T
#   E = return_gradient_type(T,first(x))
#   F = return_gradient_type(E,first(x))
#   c = zeros(F,length(x))
#   CachedArray(c)
# end

# function evaluate_hessian!(c,f::ConstantField,x)
#   evaluate_gradient!(c,f,x)
# end


# Make Function behave like Field

function_field(f::Function) =  GenericField(f)

const FunctionField{F} = GenericField{F} where F<:Function

function return_gradient_cache(f::FunctionField,x)
  gradient(f.object)
end

@inline evaluate_gradient!(c,f::FunctionField,x) = c(x)

# function return_cache(f::FunctionField,x::AbstractArray{<:Point})
#   nx = length(x)
#   Te = eltype(x)
#   c = zeros(return_type(f.f,Te),nx)
#   CachedArray(c)
# end

# function evaluate!(c,f::FunctionField,x::AbstractArray{<:Point})
#   nx = length(x)
#   setsize!(c,(nx,))
#   for i in eachindex(x)
#     c[i] = f.f(x[i])
#   end
#   c
# end

# function evaluate!(cache,a::Function,x::AbstractArray{<:Point})
#   evaluate!(cache,FunctionField(a),x)
# end

# function return_cache(f::Function,x::AbstractArray{<:Point})
#   return_cache(FunctionField(f),x)
# end

# Differentiation

struct Gradient{F} <: Field
  object::F
end

gradient(f::Field) = Gradient(f)

gradient(f::GenericField{Gradient}) = Hessian(f.object.object)

@inline evaluate!(cache,f::Gradient,x) = evaluate_gradient!(cache,f.object,x)

@inline return_cache(f::Gradient,x) = return_gradient_cache(f.object,x)

struct Hessian{F} <: Field
  object::F
end

gradient(f::Gradient) = Hessian(f.object)

@inline evaluate!(cache,f::Hessian,x) = evaluate_hessian!(cache,f.object,x)

@inline return_cache(f::Hessian,x) = return_hessian_cache(f.object,x)

@inline function gradient(f::Hessian)
  @unreachable "Default implementation of 3rt order derivatives not available"
end

# Composition

Base.:∘(f::Field,g::Field) = Operation(f)(g)

# Broadcast fields

struct BroadcastField{T<:Field} <: Field
  field::T
end

function return_cache(f::BroadcastField,x::AbstractArray{<:Point})
  T = return_type(f.field,first(x))
  s = size(x)
  ab = zeros(T,s)
  cb = CachedArray(ab)
  cf = return_cache(f.field,first(x))
  cb, cf
end

function return_cache(f::BroadcastField,x::Point)
  return_cache(f.field,x)
end

function evaluate!(c,f::BroadcastField,x::AbstractArray{<:Point})
  cb, cf = c
  sx = size(x)
  setsize!(cb,sx)
  for i in eachindex(x)
    @inbounds cb[i] = evaluate!(cf,f.field,x[i])
  end
  cb.array
end

gradient(b::BroadcastField) = BroadcastField(gradient(b.field))

# Testers

function test_field(
  f::FieldOrFieldArray,
  x::Tuple,
  v,cmp=(==);
  grad=nothing,
  hessian=nothing)

  x, = x

  @test isa(x,Union{Point,AbstractVector{<:Point}})

  w = evaluate(f,x)

  @test cmp(w,v)
  @test typeof(w) == return_type(f,x)

  cf = return_cache(f,x)
  r = evaluate!(cf,f,x)
  @test cmp(r,v)

  if x isa AbstractVector{<:Point}
     np, = size(w)
     @test length(x) == np

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
