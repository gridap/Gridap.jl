"""
    Chebyshev{kind} <: Polynomial

Type representing Chebyshev polynomials of first and second kind
where `kind` is either `:T` or `:U` for first and second kind Chebyshev polynomials respectively.
"""
struct Chebyshev{kind} <: Polynomial end

isHierarchical(::Chebyshev) = true

"""
    ChebyshevBasis{D,V,kind,K} = TensorPolynomialBasis{D,V,K,Chebyshev{kind}}

Multivariate scalar' or `Multivalue`'d Chebyshev basis, see [`TensorPolynomialBasis`](@ref)
"""
const ChebyshevBasis{D,V,kind,K} = TensorPolynomialBasis{D,V,K,Chebyshev{kind}}

"""
    ChebyshevBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector; kind=:T)
    ChebyshevBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function; kind=:T])
    ChebyshevBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function; kind=:T])

Convenience constructors of [`ChebyshevBasis`](@ref).
"""
ChebyshevBasis(args...; kind=:T) = TensorPolynomialBasis(Chebyshev{kind}, args...)

TensorPolynomialBasis{D}(::Type{Chebyshev{:U}}, args...) where D = @notimplemented "1D evaluation for second kind needed here"

QGradChebyshevBasis(args...; kind=:T)     = QGradBasis(Chebyshev{kind}, args...)
#PGradChebyshevBasis(args...; kind=:T)     = PGradBasis(Chebyshev{kind}, args...)
QCurlGradChebyshevBasis(args...; kind=:T) = QCurlGradBasis(Chebyshev{kind}, args...)
#PCurlGradChebyshevBasis(args...; kind=:T) = PCurlGradBasis(Chebyshev{kind}, args...)


# 1D evaluation implementation

function _evaluate_1d!(
  ::Type{Chebyshev{:T}},::Val{0},v::AbstractMatrix{T},x,d) where T<:Number

  @inbounds v[d,1] = one(T)
end

function _evaluate_1d!(
  ::Type{Chebyshev{:T}},::Val{K},v::AbstractMatrix{T},x,d) where {K,T<:Number}

  n = K + 1        # n > 1
  ξ = (2*x[d] - 1) # ξ ∈ [-1,1]
  ξ2 = 2*ξ

  @inbounds v[d,1] = one(T)
  @inbounds v[d,2] = ξ
  for i in 3:n
    @inbounds v[d,i] = v[d,i-1]*ξ2 - v[d,i-2]
  end
end


function _gradient_1d!(
  ::Type{Chebyshev{:T}},::Val{0},g::AbstractMatrix{T},x,d) where T<:Number

  @inbounds g[d,1] = zero(T)
end

function _gradient_1d!(
  ::Type{Chebyshev{:T}},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number}

  n = K + 1  # n>1
  z = zero(T)
  o = one(T)
  ξ = T(2*x[d] - 1)
  dξdx = T(2.0)

  unm1 = o
  un = 2*ξ
  @inbounds g[d,1] = z               # dT_0 = 0
  @inbounds g[d,2] = dξdx*o          # dT_1 = 1*U_0 = 1
  for i in 3:n
    @inbounds g[d,i] = dξdx*(i-1)*un # dT_i = i*U_{i-1}
    un, unm1 = 2*ξ*un - unm1, un
  end
end


function _hessian_1d!(
  ::Type{Chebyshev{:T}},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number}

  @notimplemented
end

############################################################################################

"""
    struct QGradChebyshevPolynomialBasis{...} <: AbstractVector{Chebyshev{:T}}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Nedelec reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QGradChebyshevPolynomialBasis{D,T} <: AbstractVector{Chebyshev{:T}}
  order::Int
  terms::CartesianIndices{D}
  perms::Matrix{Int}
  function QGradChebyshevPolynomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    new{D,T}(order,terms,perms)
  end
end

Base.size(a::QGradChebyshevPolynomialBasis) = (_ndofs_qgrad_ch(a),)
Base.getindex(a::QGradChebyshevPolynomialBasis,i::Integer) = Chebyshev{:T}()
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
    num_terms(f::QGradChebyshevPolynomialBasis{D}) where {D}
"""
num_terms(f::QGradChebyshevPolynomialBasis{D}) where {D} = length(f.terms)*D

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
    K = Val(order)
    _evaluate_1d!(Chebyshev{:T},K,c,x,d)
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
    K = Val(order)
    _derivatives_1d!(Chebyshev{:T},K,(c,g),x,d)
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
struct QCurlGradChebyshevPolynomialBasis{D,T} <: AbstractVector{Chebyshev{:T}}
  qgrad::QGradChebyshevPolynomialBasis{D,T}
  function QCurlGradChebyshevPolynomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    qgrad = QGradChebyshevPolynomialBasis(T,order,terms,perms)
    new{D,T}(qgrad)
  end
end

Base.size(a::QCurlGradChebyshevPolynomialBasis) = (length(a.qgrad),)
Base.getindex(a::QCurlGradChebyshevPolynomialBasis,i::Integer) = Chebyshev{:T}()
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

