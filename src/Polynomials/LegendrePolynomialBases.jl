"""
    LegendrePType{K} <: PolynomialType{K}

Type representing Legendre polynomials of order up to `K`
"""
struct LegendrePType{K} <: PolynomialType{K} end
isHierarchical(::LegendrePType) = true

struct LegendrePolynomial <: Field end

struct LegendrePolynomialBasis{D,T} <: AbstractVector{LegendrePolynomial}
  orders::NTuple{D,Int}
  terms::Vector{CartesianIndex{D}}
  function LegendrePolynomialBasis{D}(
    ::Type{T}, orders::NTuple{D,Int}, terms::Vector{CartesianIndex{D}}) where {D,T}
    new{D,T}(orders,terms)
  end
end

@inline Base.size(a::LegendrePolynomialBasis{D,T}) where {D,T} = (length(a.terms)*num_indep_components(T),)
@inline Base.getindex(a::LegendrePolynomialBasis,i::Integer) = LegendrePolynomial()
@inline Base.IndexStyle(::LegendrePolynomialBasis) = IndexLinear()

function LegendrePolynomialBasis{D}(
  ::Type{T}, orders::NTuple{D,Int}, filter::Function=_q_filter) where {D,T}

  terms = _define_terms(filter, orders)
  LegendrePolynomialBasis{D}(T,orders,terms)
end

function LegendrePolynomialBasis{D}(
  ::Type{T}, order::Int, filter::Function=_q_filter) where {D,T}

  orders = tfill(order,Val{D}())
  LegendrePolynomialBasis{D}(T,orders,filter)
end

# API

function get_exponents(b::LegendrePolynomialBasis)
  indexbase = 1
  [Tuple(t) .- indexbase for t in b.terms]
end

function get_order(b::LegendrePolynomialBasis)
  maximum(b.orders)
end

function get_orders(b::LegendrePolynomialBasis)
  b.orders
end

return_type(::LegendrePolynomialBasis{D,T}) where {D,T} = T

# Field implementation

function return_cache(f::LegendrePolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f)
  n = 1 + _maximum(f.orders)
  r = CachedArray(zeros(T,(np,ndof)))
  v = CachedArray(zeros(T,(ndof,)))
  c = CachedArray(zeros(eltype(T),(D,n)))
  (r, v, c)
end

function evaluate!(cache,f::LegendrePolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = length(f)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _evaluate_nd_leg!(v,xi,f.orders,f.terms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{1,LegendrePolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}

  f = fg.fa
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f)
  xi = testitem(x)
  T = gradient_type(V,xi)
  n = 1 + _maximum(f.orders)
  r = CachedArray(zeros(T,(np,ndof)))
  v = CachedArray(zeros(T,(ndof,)))
  c = CachedArray(zeros(eltype(T),(D,n)))
  g = CachedArray(zeros(eltype(T),(D,n)))
  (r, v, c, g)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{1,LegendrePolynomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = length(f)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_leg!(v,xi,f.orders,f.terms,c,g,T)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{2,LegendrePolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}

  f = fg.fa
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f)
  xi = testitem(x)
  T = gradient_type(gradient_type(V,xi),xi)
  n = 1 + _maximum(f.orders)
  r = CachedArray(zeros(T,(np,ndof)))
  v = CachedArray(zeros(T,(ndof,)))
  c = CachedArray(zeros(eltype(T),(D,n)))
  g = CachedArray(zeros(eltype(T),(D,n)))
  h = CachedArray(zeros(eltype(T),(D,n)))
  (r, v, c, g, h)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{2,LegendrePolynomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  r, v, c, g, h = cache
  np = length(x)
  ndof = length(f)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  setsize!(h,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _hessian_nd_leg!(v,xi,f.orders,f.terms,c,g,h,T)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

# Optimizing evaluation at a single point

function return_cache(f::AbstractVector{LegendrePolynomial},x::Point)
  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(cache,f::AbstractVector{LegendrePolynomial},x::Point)
  r, cf, xs = cache
  xs[1] = x
  v = evaluate!(cf,f,xs)
  ndof = size(v,2)
  setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

function return_cache(
  f::FieldGradientArray{N,<:AbstractVector{LegendrePolynomial}}, x::Point) where {N}
  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(
  cache, f::FieldGradientArray{N,<:AbstractVector{LegendrePolynomial}}, x::Point) where {N}
  r, cf, xs = cache
  xs[1] = x
  v = evaluate!(cf,f,xs)
  ndof = size(v,2)
  setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

# Helpers

function _evaluate_1d!(
  ::Type{LegendrePType{K}}, v::AbstractMatrix{T},x,d) where {K,T<:Number}

  n = K + 1
  o = one(T)
  @inbounds v[d,1] = o
  if n > 1
    ξ = ( 2*x[d] - 1 )
    for i in 2:n
      # The sqrt(2i-1) factor normalizes the basis polynomial for L2 scalar product on ξ∈[0,1], indeed:
      # ∫[0,1] Pn(2ξ-1)^2 dξ = 1/2 ∫[-1,1] Pn(t)^2 dt = 1/(2n+1)
      # C.f. Eq. (1.25) in Section 1.1.5 in Ern & Guermond book (2013).
      @inbounds v[d,i] = sqrt(2*i-1)*jacobi(ξ,i-1,0,0)
    end
  end
end

function _gradient_1d!(
  ::Type{LegendrePType{K}}, g::AbstractMatrix{T},x,d) where {K,T<:Number}

  n = K + 1
  z = zero(T)
  @inbounds g[d,1] = z
  if n > 1
    ξ = ( 2*x[d] - 1 )
    for i in 2:n
      @inbounds g[d,i] = sqrt(2*i-1)*i*jacobi(ξ,i-2,1,1)
    end
  end
end

function _hessian_1d!(
  ::Type{LegendrePType{K}}, h::AbstractMatrix{T},x,d) where {K,T<:Number}

  n = K + 1
  z = zero(T)
  @inbounds h[d,1] = z
  if n > 1
    @inbounds h[d,2] = z
    ξ = ( 2*x[d] - 1 )
    for i in 3:n
      @inbounds h[d,i] = sqrt(2*i-1)*(i*(i+1)/2)*jacobi(ξ,i-3,2,2)
    end
  end
end

function _evaluate_nd_leg!(
  v::AbstractVector{V},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    _evaluate_1d!(LegendrePType{orders[d]},c,x,d)
  end

  o = one(T)
  k = 1

  for ci in terms

    s = o
    for d in 1:dim
      @inbounds s *= c[d,ci[d]]
    end

    k = _set_value!(v,s,k)

  end

end

function _gradient_nd_leg!(
  v::AbstractVector{G},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _derivatives_1d!(LegendrePType{orders[d]},(c,g),x,d)
  end

  z = zero(Mutable(VectorValue{D,T}))
  o = one(T)
  k = 1

  for ci in terms

    s = z
    for i in eachindex(s)
      @inbounds s[i] = o
    end
    for q in 1:dim
      for d in 1:dim
        if d != q
          @inbounds s[q] *= c[d,ci[d]]
        else
          @inbounds s[q] *= g[d,ci[d]]
        end
      end
    end

    k = _set_gradient!(v,s,k,V)

  end

end

function _hessian_nd_leg!(
  v::AbstractVector{G},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _derivatives_1d!(LegendrePType{orders[d]},(c,g,h),x,d)
  end

  z = zero(Mutable(TensorValue{D,D,T}))
  o = one(T)
  k = 1

  for ci in terms

    s = z
    for i in eachindex(s)
      @inbounds s[i] = o
    end
    for r in 1:dim
      for q in 1:dim
        for d in 1:dim
          if d != q && d != r
            @inbounds s[r,q] *= c[d,ci[d]]
          elseif d == q && d ==r
            @inbounds s[r,q] *= h[d,ci[d]]
          else
            @inbounds s[r,q] *= g[d,ci[d]]
          end
        end
      end
    end

    k = _set_gradient!(v,s,k,V)

  end

end
