struct ChebyshevPolynomial <: Field end

struct ChebyshevPolynomialBasis{D,T} <: AbstractVector{ChebyshevPolynomial}
  orders::NTuple{D,Int}
  terms::Vector{CartesianIndex{D}}
  function ChebyshevPolynomialBasis{D}(
    ::Type{T}, orders::NTuple{D,Int}, terms::Vector{CartesianIndex{D}}) where {D,T}
    new{D,T}(orders,terms)
  end
end

@inline Base.size(a::ChebyshevPolynomialBasis{D,T}) where {D,T} = (length(a.terms)*num_components(T),)
@inline Base.getindex(a::ChebyshevPolynomialBasis,i::Integer) = ChebyshevPolynomial()
@inline Base.IndexStyle(::ChebyshevPolynomialBasis) = IndexLinear()

function ChebyshevPolynomialBasis{D}(
  ::Type{T}, orders::NTuple{D,Int}, filter::Function=_q_filter) where {D,T}

  terms = _define_terms(filter, orders)
  ChebyshevPolynomialBasis{D}(T,orders,terms)
end

function ChebyshevPolynomialBasis{D}(
  ::Type{T}, order::Int, filter::Function=_q_filter) where {D,T}

  orders = tfill(order,Val{D}())
  ChebyshevPolynomialBasis{D}(T,orders,filter)
end

# API

function get_exponents(b::ChebyshevPolynomialBasis)
  indexbase = 1
  [Tuple(t) .- indexbase for t in b.terms]
end

function get_order(b::ChebyshevPolynomialBasis)
  maximum(b.orders)
end

function get_orders(b::ChebyshevPolynomialBasis)
  b.orders
end

return_type(::ChebyshevPolynomialBasis{D,T}) where {D,T} = T

# Field implementation

function return_cache(f::ChebyshevPolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f.terms)*num_components(T)
  n = 1 + _maximum(f.orders)
  r = CachedArray(zeros(T,(np,ndof)))
  v = CachedArray(zeros(T,(ndof,)))
  c = CachedArray(zeros(eltype(T),(D,n)))
  (r, v, c)
end

function evaluate!(cache,f::ChebyshevPolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = length(f.terms)*num_components(T)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _evaluate_nd_ch!(v,xi,f.orders,f.terms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{1,ChebyshevPolynomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}

  f = fg.fa
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f.terms)*num_components(V)
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
  fg::FieldGradientArray{1,ChebyshevPolynomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = length(f.terms) * num_components(T)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_ch!(v,xi,f.orders,f.terms,c,g,T)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

# Optimizing evaluation at a single point

function return_cache(f::ChebyshevPolynomialBasis{D,T},x::Point) where {D,T}
  ndof = length(f.terms)*num_components(T)
  r = CachedArray(zeros(T,(ndof,)))
  xs = [x]
  cf = return_cache(f,xs)
  r, cf, xs
end

function evaluate!(cache,f::ChebyshevPolynomialBasis{D,T},x::Point) where {D,T}
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
  f::FieldGradientArray{N,ChebyshevPolynomialBasis{D,V}}, x::Point) where {N,D,V}
  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(
  cache, f::FieldGradientArray{N,ChebyshevPolynomialBasis{D,V}}, x::Point) where {N,D,V}
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

function _evaluate_1d_ch_T!(v::AbstractMatrix{T},x,order,d) where T
  n = order + 1
  o = one(T)
  @inbounds v[d,1] = o
  if n > 1
    ξ = (2*x[d] - 1) # ξ ∈ [-1,1]
    ξ2 = 2*ξ
    v[d,2] = ξ
    for i in 3:n
      @inbounds v[d,i] = v[d,i-1]*ξ2 - v[d,i-2]
    end
  end
end

function _evaluate_1d_ch_U!(v::AbstractMatrix{T},x,order,d) where T
  n = order + 1
  o = one(T)
  @inbounds v[d,1] = o
  if n > 1
    ξ = (2*x[d] - 1) # ξ ∈ [-1,1]
    ξ2 = 2*ξ
    v[d,2] = ξ2
    for i in 3:n
      @inbounds v[d,i] = v[d,i-1]*ξ2 - v[d,i-2]
    end
  end
end

function _gradient_1d_ch_T!(v::AbstractMatrix{T},x,order,d) where T
  n = order + 1
  z = zero(T)
  o = one(T)
  dξdx = T(2.0)
  @inbounds v[d,1] = z # dT_0 = 0
  if n > 1
    ξ = T(2*x[d] - 1)
    @inbounds v[d,2] = dξdx*o # dT_1 = 1*U_0 = 1
    unm1 = o
    un = 2*ξ
    for i in 3:n
      @inbounds v[d,i] = dξdx*(i-1)*un # dT_i = i*U_{i-1}
      un, unm1 = 2*ξ*un - unm1, un
    end
  end
end

function _evaluate_nd_ch!(
  v::AbstractVector{V},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    _evaluate_1d_ch_T!(c,x,orders[d],d)
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

function _gradient_nd_ch!(
  v::AbstractVector{G},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _evaluate_1d_ch_T!(c,x,orders[d],d)
    _gradient_1d_ch_T!(g,x,orders[d],d)
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

############################################################################################

"""
    struct QGradChebyshevPolynomialBasis{...} <: AbstractVector{Monomial}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Nedelec reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QGradChebyshevPolynomialBasis{D,T} <: AbstractVector{ChebyshevPolynomial}
  order::Int
  terms::CartesianIndices{D}
  perms::Matrix{Int}
  function QGradChebyshevPolynomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    new{D,T}(order,terms,perms)
  end
end

Base.size(a::QGradChebyshevPolynomialBasis) = (_ndofs_qgrad_ch(a),)
Base.getindex(a::QGradChebyshevPolynomialBasis,i::Integer) = ChebyshevPolynomial()
Base.IndexStyle(::QGradChebyshevPolynomialBasis) = IndexLinear()

"""
    QGradChebyshevPolynomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QGradChebyshevPolynomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the curl of the  functions in this basis
is in the Q space of degree `order`.
"""
function QGradChebyshevPolynomialBasis{D}(::Type{T},order::Int) where {D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _order = order + 1
  _t = tfill(_order+1,Val{D-1}())
  t = (_order,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QGradChebyshevPolynomialBasis(T,order,terms,perms)
end

"""
    num_terms(f::QGradChebyshevPolynomialBasis{D,T}) where {D,T}
"""
num_terms(f::QGradChebyshevPolynomialBasis{D,T}) where {D,T} = length(f.terms)*D

get_order(f::QGradChebyshevPolynomialBasis) = f.order

function return_cache(f::QGradChebyshevPolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad_ch(f)
  n = 1 + f.order+1
  V = VectorValue{D,T}
  r = CachedArray(zeros(V,(np,ndof)))
  v = CachedArray(zeros(V,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  (r, v, c)
end

function evaluate!(cache,f::QGradChebyshevPolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = _ndofs_qgrad_ch(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_qgrad_ch!(v,xi,f.order+1,f.terms,f.perms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{1,QGradChebyshevPolynomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad_ch(f)
  n = 1 + f.order+1
  xi = testitem(x)
  V = VectorValue{D,T}
  G = gradient_type(V,xi)
  r = CachedArray(zeros(G,(np,ndof)))
  v = CachedArray(zeros(G,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  g = CachedArray(zeros(T,(D,n)))
  (r, v, c, g)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{1,QGradChebyshevPolynomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = _ndofs_qgrad_ch(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  V = VectorValue{D,T}
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_qgrad_ch!(v,xi,f.order+1,f.terms,f.perms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

# Helpers

_ndofs_qgrad_ch(f::QGradChebyshevPolynomialBasis{D}) where D = D*(length(f.terms))

function _evaluate_nd_qgrad_ch!(
  v::AbstractVector{V},
  x,
  order,
  terms::CartesianIndices{D},
  perms::Matrix{Int},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    _evaluate_1d_ch_T!(c,x,order,d)
  end

  o = one(T)
  k = 1
  m = zero(Mutable(V))
  js = eachindex(m)
  z = zero(T)

  for ci in terms

    for j in js

      @inbounds for i in js
        m[i] = z
      end

      s = o
      @inbounds for d in 1:dim
        s *= c[d,ci[perms[d,j]]]
      end

      m[j] = s
      v[k] = m
      k += 1

    end

  end

end

function _gradient_nd_qgrad_ch!(
  v::AbstractVector{G},
  x,
  order,
  terms::CartesianIndices{D},
  perms::Matrix{Int},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _evaluate_1d_ch_T!(c,x,order,d)
    _gradient_1d_ch_T!(g,x,order,d)
  end

  z = zero(Mutable(V))
  m = zero(Mutable(G))
  js = eachindex(z)
  mjs = eachindex(m)
  o = one(T)
  zi = zero(T)
  k = 1

  for ci in terms

    for j in js

      s = z
      for i in js
        s[i] = o
      end

      for q in 1:dim
        for d in 1:dim
          if d != q
            @inbounds s[q] *= c[d,ci[perms[d,j]]]
          else
            @inbounds s[q] *= g[d,ci[perms[d,j]]]
          end
        end
      end

      @inbounds for i in mjs
        m[i] = zi
      end

      for i in js
        @inbounds m[i,j] = s[i]
      end
        @inbounds v[k] = m
      k += 1

    end

  end

end

############################################################################################

"""
    struct QCurlGradChebyshevPolynomialBasis{...} <: AbstractArray{Monomial}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Raviart-Thomas reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QCurlGradChebyshevPolynomialBasis{D,T} <: AbstractVector{ChebyshevPolynomial}
  qgrad::QGradChebyshevPolynomialBasis{D,T}
  function QCurlGradChebyshevPolynomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    qgrad = QGradChebyshevPolynomialBasis(T,order,terms,perms)
    new{D,T}(qgrad)
  end
end

Base.size(a::QCurlGradChebyshevPolynomialBasis) = (length(a.qgrad),)
Base.getindex(a::QCurlGradChebyshevPolynomialBasis,i::Integer) = ChebyshevPolynomial()
Base.IndexStyle(::QCurlGradChebyshevPolynomialBasis) = IndexLinear()

"""
    QCurlGradChebyshevPolynomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QCurlGradChebyshevPolynomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function QCurlGradChebyshevPolynomialBasis{D}(::Type{T},order::Int) where {D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _order = order+1
  _t = tfill(_order,Val{D-1}())
  t = (_order+1,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QCurlGradChebyshevPolynomialBasis(T,order,terms,perms)
end

return_type(::QCurlGradChebyshevPolynomialBasis{D,T}) where {D,T} = T

function return_cache(f::QCurlGradChebyshevPolynomialBasis,x::AbstractVector{<:Point})
  return_cache(f.qgrad,x)
end

function evaluate!(cache,f::QCurlGradChebyshevPolynomialBasis,x::AbstractVector{<:Point})
  evaluate!(cache,f.qgrad,x)
end

function return_cache(
  fg::FieldGradientArray{N,<:QCurlGradChebyshevPolynomialBasis},
  x::AbstractVector{<:Point}) where N

  f = fg.fa
  return_cache(FieldGradientArray{N}(f.qgrad),x)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{N,<:QCurlGradChebyshevPolynomialBasis},
  x::AbstractVector{<:Point}) where N

  f = fg.fa
  evaluate!(cache,FieldGradientArray{N}(f.qgrad),x)
end

"""
    num_terms(f::QCurlGradChebyshevPolynomialBasis{D,T}) where {D,T}
"""
num_terms(f::QCurlGradChebyshevPolynomialBasis{D,T}) where {D,T} = length(f.qgrad.terms)*D

get_order(f::QCurlGradChebyshevPolynomialBasis{D,T}) where {D,T} = get_order(f.qgrad)

