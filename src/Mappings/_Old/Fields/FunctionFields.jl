function_field(f::Function) =  GenericField(f)

const FunctionField{F} = GenericField{F<:Function}

function return_cache(f::FunctionField,x::AbstractArray{<:Point})
  nx = length(x)
  Te = eltype(x)
  c = zeros(return_type(f.f,Te),nx)
  CachedArray(c)
end

function evaluate!(c,f::FunctionField,x::AbstractArray{<:Point})
  nx = length(x)
  setsize!(c,(nx,))
  for i in eachindex(x)
    c[i] = f.f(x[i])
  end
  c
end

function gradient(f::FunctionField)
  FunctionField(gradient(f.f))
end

# TODO for the moment the user specifies the gradient of a function
# like f(x) = 1, g(x) = zero(x), gradient(::typeof(f)) = g
# A better default implementation with AD will replace this once merged with Gridap

gradient(f::Function) = @abstractmethod

# Make Function and numbers behave like Fields

function evaluate!(cache,a::Function,x::AbstractArray{<:Point})
  evaluate!(cache,FunctionField(a),x)
en

function return_cache(f::Function,x::AbstractArray{<:Point})
  return_cache(FunctionField(f),x)
end

# Make array of Function and Array of number behave like array of Fields

function evaluate!(cache,a::AbstractArray{<:Function},x::Point)
  map(f->evaluate(f,x),a) # TODO cache
end

function evaluate!(cache,a::AbstractArray{<:Function},x::AbstractVector{<:Point})
  fx = map(f->evaluate(f,x),a) # TODO cache
  T = eltype(map(eltype,fx))
  r = zeros(T,length(x),length(a)) # TODO avoid allocations with cache
  @inbounds for i in 1:length(x)
    for j in 1:length(a)
      r[i,j] = fx[j][i]
    end
  end
  r
end
