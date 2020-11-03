
"""
    struct Homothecy{D,T} <:Field
      origin::Point{D,T}
      scaling::Point{D,T}
    end

This type is slightly more general than a homothecy since we allow anisotropic scalings.

Implements the `Field` interface up to first order derivatives
"""
struct Homothecy{D,T} <:Field
  origin::Point{D,T}
  scaling::Point{D,T}
end

function return_cache(f::Homothecy,x)
  y = copy(x)
  CachedArray(y)
end

@inline function evaluate!(cache,f::Homothecy,x)
  setsize!(cache,size(x))
  y = cache.array
  @inbounds for i in eachindex(x)
    xi = x[i]
    yi = _apply_homothecy(f,xi)
    y[i] = yi
  end
  y
end

function _apply_homothecy(h,x)
  y = zero(mutable(x))
  t = h.origin
  s = h.scaling
  for i in 1:length(y)
    y[i] = s[i]*(x[i]-t[i]) + t[i]
  end
  Point(y)
end

struct HomothecyGrad{D,T} <: Field
  scaling::Point{D,T}
end

function field_gradient(h::Homothecy)
  HomothecyGrad(h.scaling)
end

function return_cache(f::HomothecyGrad,x)
  xi = testitem(x)
  T = typeof(xi)
  G = return_gradient_type(T,xi)
  j = similar(x,G)
  CachedArray(j)
end

function evaluate!(cache,f::HomothecyGrad,x)
  setsize!(cache,size(x))
  y = cache.array
  G = eltype(y)
  yi = _apply_homothecy_gradient(f,G)
  @inbounds for i in eachindex(y)
    y[i] = yi
  end
  y
end

function _apply_homothecy_gradient(h,G)
  y = zero(mutable(G))
  s = h.scaling
  for i in 1:length(s)
    y[i,i] = s[i]
  end
  TensorValue(y)
end
