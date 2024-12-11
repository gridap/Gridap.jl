"""
    Monomial <: Polynomial

Type representing the monomial polynomials
"""
struct Monomial <: Polynomial   end

isHierarchical(::Monomial) = true

"""
    MonomialBasis{D,T,K} = TensorPolynomialBasis{D,T,K,Monomial}

Multivariate scalar' or `Multivalue`'d monomial basis, see [`TensorPolynomialBasis`](@ref)
"""
const MonomialBasis{D,T,K} = TensorPolynomialBasis{D,T,K,Monomial}

"""
    MonomialBasis{D}(::Type{T}, order::Int, terms::Vector) where {D,T}
    MonomialBasis{D}(::Type{T}, orders::Tuple [, filter::Function]) where {D,T}
    MonomialBasis{D}(::Type{T}, order::Int [, filter::Function]) where {D,T}

Convenience constructors of MonomialBasis{D,T}.
"""
MonomialBasis{D}(args...) where {D} = TensorPolynomialBasis{D}(Monomial, args...)

# 1D evaluation implementation

function _evaluate_1d!(::Type{Monomial}, k, v::AbstractMatrix{T},x,d) where T<:Number
  n = k + 1
  @inbounds xd = x[d]
  xn = one(T)
  for i in 1:n
    @inbounds v[d,i] = xn
    xn *= xd
  end
end

function _gradient_1d!(::Type{Monomial}, k, g::AbstractMatrix{T},x,d) where T<:Number
  n = k + 1
  z = zero(T)
  @inbounds g[d,1] = z
  @inbounds xd = x[d]
  xn = one(T)
  for i in 2:n
    @inbounds g[d,i] = (i-1)*xn
    xn *= xd
  end
end

function _hessian_1d!(::Type{Monomial}, k, h::AbstractMatrix{T},x,d) where T<:Number
  n = k + 1
  z = zero(T)
  @inbounds h[d,1] = z
  if n>1
    @inbounds h[d,2] = z
  end
  @inbounds xd = x[d]
  xn = one(T)
  for i in 3:n
    @inbounds h[d,i] = (i-1)*(i-2)*xn
    xn *= xd
  end
end

