
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

