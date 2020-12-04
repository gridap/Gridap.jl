
"""
A Field with this form
y = x⋅G + y0
"""
struct AffineMap{D,T,L} <:Field
  gradient::TensorValue{D,D,T,L}
  origin::Point{D,T}
  function AffineMap(gradient::TensorValue{D,D,T,L}, origin::Point{D,T}) where {D,T,L}
    new{D,T,L}(gradient,origin)
  end
end

@inline function evaluate!(cache,f::AffineMap,x::Point)
  G = f.gradient
  y0 = f.origin
  x⋅G + y0
end

function return_cache(f::AffineMap,x::AbstractVector{<:Point})
  T = return_type(f,testitem(x))
  y = similar(x,T,size(x))
  CachedArray(y)
end

@inline function evaluate!(cache,f::AffineMap,x::AbstractVector{<:Point})
  setsize!(cache,size(x))
  y = cache.array
  G = f.gradient
  y0 = f.origin
  for i in eachindex(x)
    xi = x[i]
    yi = xi⋅G + y0
    y[i] = yi
  end
  y
end

function gradient(h::AffineMap)
  ConstantField(h.gradient)
end

