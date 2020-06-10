
"""
"""
struct AffineMap{D,T,L} <:Field
  jacobian::TensorValue{D,D,T,L}
  origin::Point{D,T}
  function AffineMap(jacobian::TensorValue{D,D,T,L}, origin::Point{D,T}) where {D,T,L}
    new{D,T,L}(jacobian,origin)
  end
end

function field_cache(f::AffineMap,x)
  y = copy(x)
  CachedArray(y)
end

@inline function evaluate_field!(cache,f::AffineMap,x)
  setsize!(cache,size(x))
  y = cache.array
  @inbounds for i in eachindex(x)
    xi = x[i]
    yi = _apply_affine_map(f,xi)
    y[i] = yi
  end
  y
end

function _apply_affine_map(h,x)
  t = h.origin
  s = h.jacobian
  (s⋅x)+t
end

struct AffineMapGrad{D,T,L} <: Field
  jacobian::TensorValue{D,D,T,L}
end

function field_gradient(h::AffineMap)
  AffineMapGrad(h.jacobian)
end

function field_cache(f::AffineMapGrad,x)
  xi = testitem(x)
  T = typeof(xi)
  G = gradient_type(T,xi)
  j = similar(x,G)
  CachedArray(j)
end

function evaluate_field!(cache,f::AffineMapGrad,x)
  setsize!(cache,size(x))
  y = cache.array
  G = eltype(y)
  yi = f.jacobian
  @inbounds for i in eachindex(y)
    y[i] = yi
  end
  y
end

