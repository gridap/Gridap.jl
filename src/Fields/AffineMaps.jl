
"""
A Field with this form
y = x⋅G + y0
"""
struct AffineMap{D1,D2,T,L} <:Field
  gradient::TensorValue{D1,D2,T,L}
  origin::Point{D2,T}
  function AffineMap(
    gradient::TensorValue{D1,D2,T,L},
    origin::Point{D2,T}) where {D1,D2,T,L}

    new{D1,D2,T,L}(gradient,origin)
  end
end

affine_map(gradient,origin) = AffineMap(gradient,origin)

function evaluate!(cache,f::AffineMap,x::Point)
  G = f.gradient
  y0 = f.origin
  x⋅G + y0
end

function return_cache(f::AffineMap,x::AbstractVector{<:Point})
  T = return_type(f,testitem(x))
  y = similar(x,T,size(x))
  CachedArray(y)
end

function evaluate!(cache,f::AffineMap,x::AbstractVector{<:Point})
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

#function push_∇∇(∇∇a::Number,Jt_inv::MultiValue{Tuple{D,D}} where D)
#  #Jt_inv⋅Jt_inv⋅∇∇a
#  Jt_inv⋅∇∇a⋅transpose(Jt_inv)
#end

function push_∇∇(∇∇a::Number,Jt_inv::MultiValue{Tuple{D1,D2}} where {D1,D2})
  _permdims_for_∇∇(Jt_inv⋅_permdims_for_∇∇(∇∇a)⋅transpose(Jt_inv))
end

function _permdims_for_∇∇(a::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  a
end
@generated function _permdims_for_∇∇(a::MultiValue{Tuple{D1,D2,D3}}) where {D1,D2,D3}
  ss = String[]
  for k in 1:D2
    for j in 1:D3
      for i in 1:D1
        push!(ss,"a[$i,$k,$j],")
      end
    end
  end
  str =  join(ss)
  Meta.parse("ThirdOrderTensorValue{$D1,$D3,$D2}($str)")
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

function inverse_map(f::AffineMap)
  Jt = f.gradient
  y0 = f.origin
  invJt = pinvJt(Jt)
  x0 = -y0⋅invJt
  AffineMap(invJt,x0)
end

function lazy_map(::typeof(∇),a::LazyArray{<:Fill{typeof(affine_map)}})
  gradients = a.args[1]
  lazy_map(constant_field,gradients)
end

function Base.zero(::Type{<:AffineMap{D1,D2,T}}) where {D1,D2,T}
  gradient = TensorValue{D1,D2}(tfill(zero(T),Val{D1*D2}()))
  origin = Point{D2,T}(tfill(zero(T),Val{D2}()))
  AffineMap(gradient,origin)
end
