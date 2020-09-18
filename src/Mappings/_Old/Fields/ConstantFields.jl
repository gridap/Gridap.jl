const ConstantField{T} = GenericField{T<:Number}

constant_field(v::Number) =  ConstantField(v)

function return_cache(f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  c = zeros(typeof(f.v),nx)
  CachedArray(c)
end

function evaluate!(c,f::ConstantField,x::AbstractArray{<:Point})
  nx = length(x)
  setsize!(c,(nx,))
  r = c.array
  for i in eachindex(x)
    @inbounds r[i] = f.v
  end
  r
end

# Gradients

gradient(f::ConstantField) = GenericField(Gradient(f.object))

function return_gradient_cache(f::ConstantField{T},x) where T
  E = return_gradient_type(T,first(x))
  c = zeros(E,length(x))
  CachedArray(c)
end

function evaluate_gradient!(c,f::ConstantField,x)
  nx = length(x)
  if size(c) != nx
    setsize!(c,(nx,))
    c .= zero(eltype(c))
  end
  c
end

function return_hessian_cache(f::ConstantField{T},x) where T
  E = return_gradient_type(T,first(x))
  F = return_gradient_type(E,first(x))
  c = zeros(F,length(x))
  CachedArray(c)
end

function evaluate_hessian!(c,f::ConstantField,x)
  evaluate_gradient!(c,f,x)
end

# Array

const ConstantFieldArray{F} = GenericField{F<:AbstractArray{<:Number}}

function return_cache(f::ConstantFieldArray,x::AbstractArray{<:Point})
  nx = length(x)
  sv = size(f)
  s = (nx,sv...)
  c = zeros(T,s)
  CachedArray(c)
end

function evaluate!(c,f::ConstantFieldArray,x::AbstractArray{<:Point})
  nx = length(x)
  sv = size(f)
  s = (nx,sv...)
  setsize!(c,s)
  cis = CartesianIndices(f)
  r = c.array
  for i in eachindex(x)
    for ci in cis
      @inbounds r[i,ci] = f[ci].v
    end
  end
  r
end

# Gradient of arrays

function return_gradient_cache(f::ConstantFieldArray,x)
  T = eltype(f.object)
  E = return_gradient_type(T,first(x))
  s = (length(x), size(f)...)
  c = zeros(E,s)
  CachedArray(c)
end

function evaluate_gradient!(c,f::ConstantFieldArray,x)
  nx = length(x)
  if size(c) != nx
    s = (nx,size(f)...)
    setsize!(c,s)
    c .= zero(eltype(c))
  end
  c
end

function return_hessian_cache(f::ConstantFieldArray,x)
  T = eltype(f.object)
  E = return_gradient_type(T,first(x))
  F = return_gradient_type(E,first(x))
  s = (length(x), size(f)...)
  c = zeros(F,s)
  CachedArray(c)
end

function evaluate_hessian!(c,f::ConstantFieldArray,x)
  nx = length(x)
  if size(c) != nx
    s = (nx,size(f)...)
    setsize!(c,s)
    c .= zero(eltype(c))
  end
  c
end

# Make Numbers behave like Fields

function evaluate!(cache,a::Number,x::AbstractArray{<:Point})
  evaluate!(cache,ConstantField(a),x)
en

function return_cache(a::Number,x::AbstractArray{<:Point})
  return_cache(ConstantField(f),x)
end

# Make Arrays behave like Fields

function evaluate!(cache,a::AbstractArray{<:Number},x::AbstractArray{<:Point})
  evaluate!(cache,ConstantFieldArray(a),x)
en

function return_cache(a::AbstractArray{<:Number},x::AbstractArray{<:Point})
  return_cache(ConstantFieldArray(f),x)
end

function evaluate_gradient!(cache,a::AbstractArray{<:Number},x::Point)
  evaluate_gradient!(cache,ConstantFieldArray(a),x::Point)
end
