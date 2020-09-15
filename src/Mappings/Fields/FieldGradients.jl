# Default gradient and hessian implementation

# """
# $(SIGNATURES)

# Returns another field that represents the gradient of the given one
# """
function gradient(f::NewField)
  FieldGradient(f)
end

function gradient(af::FieldArray)
  FieldGradientArray(af)
end

struct FieldGradient{F<:NewField} <: NewField
  field::F
  FieldGradient(f) = new{typeof(f)}(f)
end

@inline function evaluate!(cache,f::FieldGradient,x)
  evaluate_gradient!(cache,f.field,x)
end

@inline function return_cache(f::FieldGradient,x::AbstractArray{<:Point})
  return_gradient_cache(f.field,x)
end

@inline function gradient(f::FieldGradient)
  FieldHessian(f.field)
end

struct FieldHessian{F} <: NewField
  field::F
  FieldHessian(f) = new{typeof(f)}(f)
end

@inline function evaluate!(cache,f::FieldHessian,x::AbstractArray{<:Point})
  evaluate_hessian!(cache,f.field,x)
end

@inline function return_cache(f::FieldHessian,x::AbstractArray{<:Point})
  return_hessian_cache(f.field,x)
end

@inline function gradient(f::FieldHessian)
  @unreachable "Default implementation of 3rt order derivatives not available"
end

# Array of gradients

function gradient(f::FieldArray)
  FieldGradientArray(f)
end

struct FieldGradientArray{F<:NewField,N,A<:AbstractArray{F,N}} <: AbstractArray{FieldGradient{F},N}
  i_to_f::A
end

@inline function evaluate!(cache,f::FieldGradientArray,x::AbstractArray{<:Point})
  evaluate_gradient!(cache,f.i_to_f,x)
end

@inline function return_cache(f::FieldGradientArray,x::AbstractArray{<:Point})
  return_gradient_cache(f.i_to_f,x)
end

Base.size(a::FieldGradientArray) = size(a.i_to_f)
Base.getindex(a::FieldGradientArray{F,N},i::Vararg{Int,N}) where {F,N} = FieldGradient(a.i_to_f[i...])
Base.IndexStyle(::Type{<:FieldGradientArray{F,N,A}}) where {F,N,A} = IndexStyle(A)

function gradient(f::FieldGradientArray)
  FieldHessianArray(f.i_to_f)
end

struct FieldHessianArray{F<:NewField,N,A<:AbstractArray{F,N}} <: AbstractArray{FieldHessian{F},N}
  i_to_f::A
end

@inline function evaluate!(cache,f::FieldHessianArray,x::AbstractArray{<:Point})
  evaluate_hessian!(cache,f.i_to_f,x)
end

@inline function return_cache(f::FieldHessianArray,x::AbstractArray{<:Point})
  return_hessian_cache(f.i_to_f,x)
end

Base.size(a::FieldHessianArray) = size(a.i_to_f)
Base.getindex(a::FieldHessianArray{F,N},i::Vararg{Int,N}) where {F,N} = FieldHessian(a.i_to_f[i...])
Base.IndexStyle(::Type{<:FieldHessianArray{F,N,A}}) where {F,N,A} = IndexStyle(A)

# General implementation not optimized

# function return_cache(i_to_f::FieldArrayGradient,x::AbstractVector{<:Point})
#   i_to_fx = [evaluate(f,x) for f in i_to_f]
#   i_to_T = map( eltype , i_to_fx )
#   T = eltype([ zero(Ti) for Ti in i_to_T ]) # OPTION 1
#   T = promote_type(i_to_T...) # OPTION 2
#   r = zeros(T,(length(x),size(i_to_f)...))
#   cr = CachedArray(r)
#   i_to_cache = [return_cache(f,x) for t in i_to_f]
#   cr, i_to_cache
# end

# function evaluate!(cache,i_to_f::FieldArrayGradient,x::AbstractVector{<:Point})
#   cr, i_to_cache = cache
#   setsize!(cr,(length(x),size(i_to_f)...))
#   r = cr.array
#   for i in eachindex(i_to_f)
#    fx = evaluate!(i_to_cache[i],i_to_f[i],x)
#    for p in 1:length(x)
#      r[p,i] = fx[p]
#    end
#   end
#   r
# end
