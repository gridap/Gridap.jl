const Point{D,T} = VectorValue{D,T}

abstract type NewField <: Mapping end

const FieldVector{T} = AbstractVector{T} where T<:NewField

const FieldArray{T,N} = AbstractArray{T,N} where {T<:NewField,N}

const FieldOrFieldArray = Union{NewField,FieldArray}

function evaluate!(cache,f::NewField,x::Point)
  @notimplemented
end

function evaluate!(cache,f::NewField,x::AbstractArray{<:Point})
  @notimplemented
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
