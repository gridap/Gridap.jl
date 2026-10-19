
"""
    struct AffineMap <: Map
"""
struct AffineMap <: Map end

function evaluate!(cache,::AffineMap,G::TensorValue{D1,D2},y0::Point{D2},x::Point{D1}) where {D1,D2}
  x⋅G + y0
end

"""
    struct AffineField{D1,D2,T,L} <: Field

A `Field` with the form:

    y = x⋅G + y0

with `G::TensorValue{D1,D2,T,L}` and `y0::Point{D2,T}`.
"""
struct AffineField{D1,D2,T,L} <: Field
  gradient::TensorValue{D1,D2,T,L}
  origin::Point{D2,T}
  function AffineField(
    gradient::TensorValue{D1,D2,T,L},
    origin::Point{D2,T}
  ) where {D1,D2,T,L}
    new{D1,D2,T,L}(gradient,origin)
  end
end

"""
    affine_map(gradient, origin) = AffineField(gradient, origin)

See [`AffineField`](@ref).
"""
affine_map(gradient,origin) = AffineField(gradient,origin)

function Base.zero(::Type{<:AffineField{D1,D2,T}}) where {D1,D2,T}
  gradient = TensorValue{D1,D2}(tfill(zero(T),Val{D1*D2}()))
  origin = Point{D2,T}(tfill(zero(T),Val{D2}()))
  AffineField(gradient,origin)
end

function evaluate!(cache,f::AffineField,x::Point)
  G = f.gradient
  y0 = f.origin
  x⋅G + y0
end

function return_cache(f::AffineField,x::AbstractVector{<:Point})
  T = return_type(f,testitem(x))
  y = similar(x,T,size(x))
  CachedArray(y)
end

function evaluate!(cache,f::AffineField,x::AbstractVector{<:Point})
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

function gradient(h::AffineField)
  ConstantField(h.gradient)
end

function push_∇∇(∇∇a::Field,ϕ::AffineField)
  # Assuming ϕ is affine map
  Jt = ∇(ϕ)
  Jt_inv = pinvJt(Jt)
  Operation(push_∇∇)(∇∇a, Jt_inv)
end

function push_∇∇(∇∇a::Number,Jt_inv::MultiValue{Tuple{D1,D2}} where {D1,D2})
  push_∇ⁿ(∇∇a, Jt_inv, Val(2))
end

# Pushforward of the N-th gradient by an affine map
#
# For x = Jξ + b, ∂ₓ = J⁻ᵀ ∂_ξ in every derivative index, so the physical N-th gradient is the
# reference one with each of its first N indices contracted with J⁻ᵀ (`pinvJt`).

function push_∇ⁿ(∇ⁿa::Field,ϕ::AffineField,::Val{N}) where N
  Jt = ∇(ϕ)
  Jt_inv = pinvJt(Jt)
  Operation(_PushNthGradientMap{N}())(∇ⁿa, Jt_inv)
end

"""
    push_∇ⁿ(∇ⁿa::MultiValue, Jt_inv::MultiValue, Val(N))

Pushforward of the reference `N`-th gradient `∇ⁿa` by an affine map with inverse transposed
Jacobian `Jt_inv`: each of the first `N` indices is contracted with `Jt_inv`.
"""
function push_∇ⁿ(∇ⁿa::MultiValue,Jt_inv::MultiValue,::Val{N}) where N
  # Contracting the last derivative index moves the new index to the front,
  # so after N rounds the derivative indices are back in their original order.
  t = ∇ⁿa
  for _ in 1:N
    t = tensor_contraction(Jt_inv, t, Val(2), Val(N))
  end
  t
end

struct _PushNthGradientMap{N} <: Function end
(::_PushNthGradientMap{N})(t, Jt_inv) where N = push_∇ⁿ(t, Jt_inv, Val(N))

function lazy_map(
  k::Broadcasting{typeof(push_∇∇)},
  cell_∇∇a::AbstractArray,
  cell_map::AbstractArray{<:AffineField})
  cell_Jt = lazy_map(∇,cell_map)
  cell_invJt = lazy_map(Operation(pinvJt),cell_Jt)
  return lazy_map(Broadcasting(Operation(push_∇∇)),cell_∇∇a,cell_invJt)
end

function lazy_map(
  k::Broadcasting{typeof(push_∇∇)},
  cell_∇∇at::LazyArray{<:Fill{typeof(transpose)}},
  cell_map::AbstractArray{<:AffineField})
  cell_∇∇a = cell_∇∇at.args[1]
  cell_∇∇b = lazy_map(k,cell_∇∇a,cell_map)
  return lazy_map(transpose,cell_∇∇b)
end

function lazy_map(
  k::Broadcasting{PushNthGradient{N}},
  cell_∇ⁿa::AbstractArray,
  cell_map::AbstractArray{<:AffineField}) where N
  cell_Jt = lazy_map(∇,cell_map)
  cell_invJt = lazy_map(Operation(pinvJt),cell_Jt)
  return lazy_map(Broadcasting(Operation(_PushNthGradientMap{N}())),cell_∇ⁿa,cell_invJt)
end

function lazy_map(
  k::Broadcasting{PushNthGradient{N}},
  cell_∇ⁿat::LazyArray{<:Fill{typeof(transpose)}},
  cell_map::AbstractArray{<:AffineField}) where N
  cell_∇ⁿa = cell_∇ⁿat.args[1]
  return lazy_map(transpose, lazy_map(k,cell_∇ⁿa,cell_map))
end

function inverse_map(f::AffineField)
  Jt = f.gradient
  y0 = f.origin
  invJt = pinvJt(Jt)
  x0 = -y0⋅invJt
  AffineField(invJt,x0)
end

function lazy_map(::typeof(∇),a::LazyArray{<:Fill{typeof(affine_map)}})
  gradients = a.args[1]
  lazy_map(constant_field,gradients)
end

function lazy_map(
  ::typeof(evaluate),a::LazyArray{<:Fill{typeof(affine_map)}},x::AbstractArray
)
  gradients = a.args[1]
  origins = a.args[2]
  lazy_map(Broadcasting(AffineMap()),gradients,origins,x)
end

# Constructor from a simplex given by D1+1 points
function affine_map(points::NTuple{D,Point{D2,T}}) where {D,D2,T}
  origin, pk = first_and_tail(points)
  gradient = TensorValues.tensor_from_rows(map(p -> p-origin, pk))
  AffineField(gradient,origin)
end
