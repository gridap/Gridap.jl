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

# Integrate

struct Integrate <: Mapping end

function return_cache(k::Integrate,f::AbstractVector,w,j)
  T = _integrate_rt(f,w,j)
  zero(T)
end

@noinline function evaluate!(z,k::Integrate,f::AbstractVector,w,j)
  _integrate_checks(f,w,j)
  r = z
  for p in eachindex(f)
    @inbounds r = r + f[p]*w[p]*meas(j[p])
  end
  r
end

function return_cache(k::Integrate,f::AbstractArray,w,j)
  T = _integrate_rt(f,w,j)
  r = zeros(T,size(f)[2:end])
  c = CachedArray(r)
end

@inline function evaluate!(c,k::Integrate,f::AbstractArray,w,j)
  _integrate_checks(f,w,j)
  _s = size(f)
  np = _s[1]
  s = _s[2:end]
  cis = CartesianIndices(s)
  setsize!(c,s)
  z = zero(eltype(c))
  r = c.array
  for i in cis
    @inbounds r[i] = z
  end
  for p in 1:np
    @inbounds dV = w[p]*meas(j[p])
    for i in cis
      @inbounds r[i] += f[p,i]*dV
    end
  end
  r
end

function _integrate_rt(f,w,j)
  Tf = eltype(f)
  Tw = eltype(w)
  Tj = eltype(j)
  return_type(*,Tf,Tw,return_type(meas,Tj))
end

function _integrate_checks(f,w,j)
  @assert _integrate_valid_sizes(f,w,j) "integrate: sizes  mismatch."
end

function _integrate_valid_sizes(f,w,j)
  nf, = size(f)
  nw = length(w)
  nj = length(j)
  (nf == nw) && (nw == nj)
end
