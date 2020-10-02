struct MatMul <:Mapping end

@inline function return_cache(::MatMul,a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  CachedArray(a*b)
end

# @santiagobadia : Probably a more clever way to implement this...
@inline function evaluate!(c,::MatMul,a::AbstractMatrix,b::AbstractMatrix)
  setsize!(c,(size(a,1),size(b,2)))
  _inplace_mul!(c.array,a,b)
  c.array
end

@inline function evaluate!(c,::MatMul,a::AbstractMatrix,b::AbstractVector)
  setsize!(c,(size(a,1),))
  _inplace_mul!(c.array,a,b)
  c.array
end

# idem broadcast, just for completion
@inline function evaluate!(c,::MatMul,a::AbstractVector,b::AbstractMatrix)
  setsize!(c,(size(a,1),size(b,2)))
  _inplace_mul!(c.array,a,b)
  c.array
end

# @santiagobadia : I don't understand why allocations in mul! with VectorValue
# So, I have created mine, which involves 0 allocations
@inline function _inplace_mul!(r,a::AbstractMatrix,b::AbstractMatrix)
  for i in 1:size(a,1)
    for j in 1:size(b,2)
      @inbounds r[i,j] = view(a,i,:)⋅view(b,:,j)
    end
  end
end

@inline function _inplace_mul!(r,a::AbstractMatrix,b::AbstractVector)
  for i in 1:size(a,1)
    @inbounds r[i] = view(a,i,:)⋅b
  end
end

# @santiagobadia : Not so happy with this... Probably a new type of array that
# takes into account the insertion of the first axis...
struct LinCombVal <: Mapping end

@inline function return_cache(::LinCombVal,a::AbstractMatrix,b::AbstractMatrix)
  return_cache(MatMul(), a, b)
end

@inline function return_cache(::LinCombVal,a::AbstractMatrix,b::AbstractVector)
  return_cache(MatMul(), a, b)
end

@inline function return_cache(::LinCombVal,a::AbstractVector,b::AbstractMatrix)
  return_cache(MatMul(),transpose(b),a)
end

@inline function return_cache(::LinCombVal,a::AbstractVector,b::AbstractVector)
  a⋅b
end

@inline function evaluate!(c,::LinCombVal,a::AbstractMatrix,b::AbstractMatrix)
  evaluate!(c,MatMul(),a,b)
end

@inline function evaluate!(c,::LinCombVal,a::AbstractMatrix,b::AbstractVector)
  evaluate!(c,MatMul(),a,b)
end

@inline function evaluate!(c,::LinCombVal,a::AbstractVector,b::AbstractMatrix)
  evaluate!(c,MatMul(),transpose(b),a)
end

@inline function evaluate!(c,::LinCombVal,a::AbstractVector,b::AbstractVector)
  a⋅b
end
