const Point{D,T} = VectorValue{D,T}

abstract type NewField <: Mapping end

const FieldVector{T} = AbstractVector{T} where T<:NewField

const FieldArray{T,N} = AbstractArray{T,N} where {T<:NewField,N}

const FieldOrFieldArray = Union{NewField,FieldArray}

function evaluate!(cache,f::NewField,x)
  @abstractmethod
end

# Differentiation

function gradient end

const ∇ = gradient

function evaluate_gradient!(cache,f,x)
  @abstractmethod
end

function return_gradient_cache(f,x)
  @abstractmethod
end

function evaluate_hessian!(cache,f,x)
  @abstractmethod
end

function return_hessian_cache(f,x)
  @abstractmethod
end

function return_gradient_type(::Type{T},x::Point) where T
typeof(outer(zero(x),zero(T)))
end

# Make Field behave like a collection
Base.length(::NewField) = 1
Base.size(::NewField) = ()
Base.axes(::NewField) = ()
Base.IteratorSize(::Type{<:NewField}) = Base.HasShape{0}()
Base.eltype(::Type{T}) where T<:NewField = T
Base.iterate(a::NewField) = (a,nothing)
Base.iterate(a::NewField,::Nothing) = nothing

# Make Fields behave like numbers (operations are defined below)
Base.zero(a::NewField) = ZeroField(a)

struct ZeroField{F} <: NewField
  field::F
end

evaluate!(cache,z::ZeroField,x::Point) = zero(return_type(z.field,x))

function evaluate!(cache,z::ZeroField,x::AbstractArray{<:Point})
  T = return_type(z.f,testitem(x))
  zeros(T,length(x)) # TODO cache
end

# @fverdugo : TODO Implementinf this
# Base.zero(a::Type{T}) where T<:Field
# requires return_type be implemented for return_type(::Type{F},args...)
# with is challenging in general. E.g. there is no general way of doing this for
# Function (afaik). No big deal since Base.zero(a::Type{T}) is not needed in practice.


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
