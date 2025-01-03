"""
    Monomial <: Polynomial

Type representing the monomial polynomials.
"""
struct Monomial <: Polynomial   end

isHierarchical(::Type{Monomial}) = true

"""
    MonomialBasis{D,V,K} = UniformPolyBasis{D,V,K,Monomial}

Alias for monomial Multivariate scalar' or `Multivalue`'d basis, see [`UniformPolyBasis`](@ref).
"""
const MonomialBasis{D,V,K} = UniformPolyBasis{D,V,K,Monomial}

"""
    MonomialBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector)
    MonomialBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function])
    MonomialBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function])

High level constructors of [`MonomialBasis`](@ref).
"""
MonomialBasis(args...) = UniformPolyBasis(Monomial, args...)

function PGradBasis(::Type{Monomial},::Val{D},::Type{T},order::Int) where {D,T}
  NedelecPolyBasisOnSimplex{D}(Monomial,T,order)
end

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

