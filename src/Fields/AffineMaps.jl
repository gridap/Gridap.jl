
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

function push_∇∇(∇∇a::Field,ϕ::AffineMap)
  # Assuming ϕ is affine map
  Jt = ∇(ϕ)
  Jt_inv = pinvJt(Jt)
  Operation(push_∇∇)(∇∇a, Jt_inv)
end

function push_∇∇(∇∇a::Number,Jt_inv::MultiValue{Tuple{D,D}} where D)
  Jt_inv⋅Jt_inv⋅∇∇a
end

function lazy_map(
  k::Broadcasting{typeof(push_∇∇)},
  cell_∇∇a::AbstractArray,
  cell_map::AbstractArray{<:AffineMap})
  cell_Jt = lazy_map(∇,cell_map)
  cell_invJt = lazy_map(Operation(pinvJt),cell_Jt)
  lazy_map(Broadcasting(Operation(push_∇∇)),cell_∇∇a,cell_invJt)
end

function lazy_map(
  k::Broadcasting{typeof(push_∇∇)},
  cell_∇∇at::LazyArray{<:Fill{typeof(transpose)}},
  cell_map::AbstractArray{<:AffineMap})
  cell_∇∇a = cell_∇∇at.args[1]
  cell_∇∇b = lazy_map(k,cell_∇∇a,cell_map)
  cell_∇∇bt = lazy_map(transpose,cell_∇∇b)
  cell_∇∇bt
end

