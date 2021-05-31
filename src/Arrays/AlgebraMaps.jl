
function return_cache(::typeof(*),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  c = a*b
  CachedArray(c)
end

@inline function evaluate!(cache,::typeof(*),a::AbstractMatrix{<:Number},b::AbstractVector{<:Number})
  m = size(a,1)
  setsize!(cache,(m,))
  c = cache.array
  mul!(c,a,b)
  c
end

@inline function evaluate!(cache,::typeof(*),a::AbstractMatrix{<:Number},b::AbstractMatrix{<:Number})
  m = size(a,1)
  n = size(b,2)
  setsize!(cache,(m,n))
  c = cache.array
  mul!(c,a,b)
  c
end

struct MulAddMap{T} <: Map
  α::T
  β::T
end

function return_cache(k::MulAddMap,a,b,c)
  d = a*b+c
  CachedArray(d)
end

@inline function evaluate!(cache,k::MulAddMap,a,b,c)
  setsize!(cache,size(c))
  d = cache.array
  copyto!(d,c)
  mul!(d,a,b,k.α,k.β)
  d
end

struct AddEntriesMap{F} <: Map
  combine::F
end

function evaluate!(cache,k::AddEntriesMap,A,v,i,j)
  add_entries!(k.combine,A,v,i,j)
end

function evaluate!(cache,k::AddEntriesMap,A,v,i)
  add_entries!(k.combine,A,v,i)
end

struct TouchEntriesMap <: Map end

function evaluate!(cache,k::TouchEntriesMap,A,v,i,j)
  add_entries!(+,A,nothing,i,j)
end

function evaluate!(cache,k::TouchEntriesMap,A,v,i)
  add_entries!(+,A,nothing,i)
end





