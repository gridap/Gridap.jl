struct MulMapping <: Mapping end

function return_cache(k::MulMapping,a,b)
  c = a*b
  CachedArray(c)
end

@inline function evaluate!(cache,k::MulMapping,a::AbstractMatrix,b::AbstractVector)
  m = axes(a,1)
  setaxes!(cache,(m,))
  c = cache.array
  mul!(c,a,b)
  c
end

@inline function evaluate!(cache,k::MulMapping,a::AbstractMatrix,b::AbstractMatrix)
  m = axes(a,1)
  n = axes(b,2)
  setaxes!(cache,(m,n))
  c = cache.array
  mul!(c,a,b)
  c
end

struct MulAddMapping{T} <: Mapping
  α::T
  β::T
end

function return_cache(k::MulAddMapping,a,b,c)
  d = a*b+c
  CachedArray(d)
end

@inline function evaluate!(cache,k::MulAddMapping,a,b,c)
  setaxes!(cache,axes(c))
  d = cache.array
  copyto!(d,c)
  mul!(d,a,b,k.α,k.β)
  d
end
