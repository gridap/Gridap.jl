"""
    Monomial <: Polynomial

Type representing the monomial polynomials, c.f. [Monomials](@ref) section.
"""
struct Monomial <: Polynomial   end

isHierarchical(::Type{Monomial}) = true

"""
    MonomialBasis{D,V} = CartProdPolyBasis{D,V,Monomial}

Alias for cartesian product monomial basis, scalar valued or multi-valued.
"""
const MonomialBasis{D,V} = CartProdPolyBasis{D,V,Monomial}

"""
    MonomialBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector)
    MonomialBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function])
    MonomialBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function])

High level constructors of [`MonomialBasis`](@ref).
"""
MonomialBasis(args...) = CartProdPolyBasis(Monomial, args...)

# 1D evaluation implementation

function _evaluate_1d!(::Type{Monomial},::Val{K},c::AbstractMatrix{T},x,d) where {K,T<:Number}
  n = K + 1
  xn = one(T)
  @inbounds xd = x[d]

  for i in 1:n
    @inbounds c[d,i] = xn
    xn *= xd
  end
end


function _gradient_1d!(::Type{Monomial},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number}
  n = K + 1
  z = zero(T)
  xn = one(T)
  @inbounds xd = x[d]

  @inbounds g[d,1] = z
  for i in 2:n
    @inbounds g[d,i] = (i-1)*xn
    xn *= xd
  end
end


function _hessian_1d!(::Type{Monomial},::Val{0},h::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds h[d,1] = zero(T)
end

function _hessian_1d!(::Type{Monomial},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number}
  n = K + 1 # n>1
  z = zero(T)
  xn = one(T)
  @inbounds xd = x[d]

  @inbounds h[d,1] = z
  @inbounds h[d,2] = z
  for i in 3:n
    @inbounds h[d,i] = (i-1)*(i-2)*xn
    xn *= xd
  end
end


# Optimizations for 0 to 1/2 derivatives at once

function _derivatives_1d!(::Type{Monomial},v::Val_01,t::NTuple{2},x,d)
  @inline _evaluate_1d!(Monomial, v, t[1], x, d)
  @inline _gradient_1d!(Monomial, v, t[2], x, d)
end

function _derivatives_1d!(::Type{Monomial},::Val{K},t::NTuple{2},x,d) where K
  @inbounds begin
    n = K + 1 # n > 2
    v, g = t
    T = eltype(v)

    z = zero(T)
    xn = one(T)
    xd = x[d]

    v[d,1] = xn
    g[d,1] = z
    for i in 2:n
      g[d,i] = (i-1)*xn
      xn *= xd
      v[d,i] = xn
    end
  end
end


function _derivatives_1d!(::Type{Monomial},v::Val_012,t::NTuple{3},x,d)
  @inline _evaluate_1d!(Monomial, v, t[1], x, d)
  @inline _gradient_1d!(Monomial, v, t[2], x, d)
  @inline _hessian_1d!( Monomial, v, t[3], x, d)
end

function _derivatives_1d!(::Type{Monomial},::Val{K},t::NTuple{3},x,d) where K
  @inbounds begin
    n = K + 1 # n > 2
    v, g, h = t
    T = eltype(v)

    z = zero(T)
    o = one(T)
    xd = x[d]

    v[d,1] = o;  g[d,1] = z; h[d,1] = z
    v[d,2] = xd; g[d,2] = o; h[d,2] = z

    xn  = xd
    xnn = o
    for i in 3:n
      h[d,i] = (i-1)*xnn
      xnn  = (i-1)*xn
      g[d,i] = xnn
      xn *= xd
      v[d,i] = xn
    end
  end
end

