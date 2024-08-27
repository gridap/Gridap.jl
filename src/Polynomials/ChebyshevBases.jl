
struct ChebyshevPolynomial{S,T} <: Field 
  k :: Int
end

function Fields.evaluate!(cache, f::ChebyshevPolynomial{:T,T}, x::Point{1}) where T
  ξ = T(2*x.data[1] - 1) # [-1,1] -> [0,1]
  o = one(T)
  p = 2*ξ
  f.k == 0 && return o
  f.k == 1 && return ξ

  t_n, t_nm1 = ξ, o
  for n in 2:f.k
    t_n, t_nm1 = p*t_n - t_nm1, t_n
  end
  return t_n
end

function Fields.evaluate!(cache, f::ChebyshevPolynomial{:U,T}, x::Point{1}) where T
  ξ = T(2*x.data[1] - 1) # [-1,1] -> [0,1]
  o = one(T)
  p = 2*ξ
  f.k == 0 && return o
  f.k == 1 && return p

  u_n, u_nm1 = p, o
  for n in 2:f.k
    u_n, u_nm1 = p*u_n - u_nm1, u_n
  end
  return u_n
end

Fields.return_type(::ChebyshevPolynomial{S,T}) where {S,T} = T

struct ChebyshevPolynomialBasis{S,T} <: AbstractVector{ChebyshevPolynomial{S,T}}
  order :: Int
  function ChebyshevPolynomialBasis{S}(::Type{T},order::Integer) where T
    @check S ∈ [:T,:U] "ChebyshevPolynomialBasis can be of type :T (first kind) or :U (second kind)"
    @check T <: Real "ChebyshevPolynomialBasis only supports Real types"
    new{S,T}(order)
  end
end

function ChebyshevPolynomialBasis(::Type{T},order::Integer) where T
  ChebyshevPolynomialBasis{:U}(T,order)
end

@inline Base.size(a::ChebyshevPolynomialBasis) = (a.order,)
@inline Base.IndexStyle(::ChebyshevPolynomialBasis) = IndexLinear()
Base.getindex(::ChebyshevPolynomialBasis{S,T},i::Integer) where {S,T} = ChebyshevPolynomial{S,T}(i-1)

get_exponents(b::ChebyshevPolynomialBasis) = [(t-1,) for t in 1:b.order]
get_order(b::ChebyshevPolynomialBasis) = b.order
get_orders(b::ChebyshevPolynomialBasis) = (b.order,)
Fields.return_type(::ChebyshevPolynomialBasis{S,T}) where {S,T} = T

function Fields.return_cache(f::ChebyshevPolynomialBasis{S,T}, x::AbstractVector{<:Point{1}}) where {S,T}
  np = length(x)
  ndof = f.order+1
  r = CachedArray(zeros(T,(np,ndof)))
  c = CachedArray(zeros(T,(np,)))
  return r, c
end

function Fields.evaluate!(cache, f::ChebyshevPolynomialBasis{:T,T}, x::AbstractVector{<:Point{1}}) where T
  _r, _c = cache

  np = length(x)
  ndof = f.order+1
  setsize!(_r,(np,ndof))
  setsize!(_c,(ndof,))
  r = _r.array
  c = _c.array

  @inbounds r[:,1] .= one(T)
  if f.order > 0
    for i in 1:np
      ξ = T(2*x[i][1] - 1)
      @inbounds r[i,2] = ξ
      @inbounds c[i] = 2*ξ
    end
    for n in 2:f.order
      @inbounds r[:,n+1] .= c .* r[:,n] - r[:,n-1]
    end
  end

  return r
end

function Fields.evaluate!(cache, f::ChebyshevPolynomialBasis{:U,T}, x::AbstractVector{<:Point{1}}) where T
  _r, _c = cache

  np = length(x)
  ndof = f.order+1
  setsize!(_r,(np,ndof))
  setsize!(_c,(np,))
  r = _r.array
  c = _c.array

  @inbounds r[:,1] .= one(T)
  if f.order > 0
    for i in 1:np
      ξ = T(2*x[i][1] - 1)
      @inbounds r[i,2] = 2*ξ
      @inbounds c[i] = 2*ξ
    end
    for n in 2:f.order
      @inbounds r[:,n+1] .= c .* r[:,n] .- r[:,n-1]
    end
  end

  return r
end

function return_cache(
  fg::FieldGradientArray{1,ChebyshevPolynomialBasis{:T,T}},
  x::AbstractVector{<:Point{1}}
) where T
  f = fg.fa
  np = length(x)
  ndof = f.order+1

  g = ChebyshevPolynomialBasis{:U,T}(f.order-1)
  g_cache = return_cache(g,x)

  V = gradient_type(T,testitem(x))
  r = CachedArray(zeros(V,(np,ndof)))
  return (r, g, g_cache)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{1,ChebyshevPolynomialBasis{:T,T}},
  x::AbstractVector{<:Point{1}}
) where T
  f = fg.fa
  _r, g, g_cache = cache

  # TODO : Multiply by the jacobian?
  np = length(x)
  ndof = f.order+1
  setsize!(_r,(np,ndof))
  r = r.array

  # dT_n = n*U_{n-1}
  V = eltype(r)
  r[:,1] .= zero(T)
  if f.order > 0
    gx = evaluate!(g_cache, g, x)
    for n in 1:f.order
      for i in 1:np
        @inbounds r[i,n+1] = V(n * gx[i,n])
      end
    end
  end

  return r
end

function return_cache(
  fg::FieldGradientArray{1,ChebyshevPolynomialBasis{:U,T}},
  x::AbstractVector{<:Point{1}}
) where T
  f = fg.fa
  np = length(x)
  ndof = f.order+1

  g = ChebyshevPolynomialBasis{:T,T}(f.order+1)
  g_cache = return_cache(g,x)
  f_cache = return_cache(f,x)

  V = gradient_type(T,testitem(x))
  r = CachedArray(zeros(V,(np,ndof)))
  return (r, g, g_cache, f_cache)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{1,ChebyshevPolynomialBasis{:U,T}},
  x::AbstractVector{<:Point{1}}
) where T
  f = fg.fa
  _r, g, g_cache, f_cache = cache

  # TODO : Multiply by the jacobian?
  np = length(x)
  ndof = f.order+1
  setsize!(_r,(np,ndof))
  r = _r.array

  # dU_n = ((n+1)*T_{n+1} - x * U_{n}) / (x^2 - 1)
  V = eltype(r)
  r[:,1] .= zero(T)
  if f.order > 0
    fx = evaluate!(f_cache, f, x)
    gx = evaluate!(g_cache, g, x)
    for i in 1:np
      ξ = T(2*x[i][1] - 1)
      for n in 1:f.order
        @inbounds r[i,n+1] = V(((n+1) * gx[i,n+1] - ξ * fx[i,n]) / (ξ*ξ - 1))
      end
    end
  end

  return r
end
