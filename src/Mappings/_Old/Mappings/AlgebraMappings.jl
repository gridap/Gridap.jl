struct MatVecMapping <: Mapping end

function return_cache(k::MatVecMapping,a,b)
  c = a*b
  CachedArray(c)
end

@inline function evaluate!(cache,k::MatVecMapping,a::AbstractMatrix,b::AbstractVector)
  m = axes(a,1)
  setaxes!(cache,(m,))
  c = cache.array
  mul!(c,a,b)
  c
end

@inline function evaluate!(cache,k::MatVecMapping,a::AbstractMatrix,b::AbstractMatrix)
  m = axes(a,1)
  n = axes(b,2)
  setaxes!(cache,(m,n))
  c = cache.array
  mul!(c,a,b)
  c
end
