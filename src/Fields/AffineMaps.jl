
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

function gradient(h::AffineMap)
  ConstantField(h.gradient)
end

#function return_cache(f::AffineMap,x::AbstractArray{<:Point})
#  y = copy(x)
#  CachedArray(y)
#end
#
#@inline function evaluate!(cache,f::AffineMap,x)
#  setsize!(cache,size(x))
#  y = cache.array
#  @inbounds for i in eachindex(x)
#    xi = x[i]
#    yi = _apply_affine_map(f,xi)
#    y[i] = yi
#  end
#  y
#end
#
#function _apply_affine_map(h,x)
#  t = h.origin
#  s = h.jacobian
#  (s⋅x)+t
#end

#struct AffineMapGrad{D,T,L} <: Field
#  jacobian::TensorValue{D,D,T,L}
#end
#
#function gradient(h::AffineMap)
#  AffineMapGrad(h.jacobian)
#end
#
#function return_cache(f::AffineMapGrad,x)
#  xi = testitem(x)
#  T = typeof(xi)
#  G = return_gradient_type(T,xi)
#  j = similar(x,G)
#  CachedArray(j)
#end
#
#function evaluate!(cache,f::AffineMapGrad,x)
#  setsize!(cache,size(x))
#  y = cache.array
#  G = eltype(y)
#  yi = f.jacobian
#  @inbounds for i in eachindex(y)
#    y[i] = yi
#  end
#  y
#end
