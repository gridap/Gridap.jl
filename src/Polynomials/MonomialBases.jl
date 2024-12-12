"""
    Monomial <: Polynomial

Type representing the monomial polynomials
"""
struct Monomial <: Polynomial   end

isHierarchical(::Monomial) = true

"""
    MonomialBasis{D,V,K} = TensorPolynomialBasis{D,V,K,Monomial}

Multivariate scalar' or `Multivalue`'d monomial basis, see [`TensorPolynomialBasis`](@ref)
"""
const MonomialBasis{D,V,K} = TensorPolynomialBasis{D,V,K,Monomial}

"""
    MonomialBasis{D}(::Type{V}, order::Int, terms::Vector) where {D,V}
    MonomialBasis{D}(::Type{V}, orders::Tuple [, filter::Function]) where {D,V}
    MonomialBasis{D}(::Type{V}, order::Int [, filter::Function]) where {D,V}

Convenience constructors of MonomialBasis{D,V}.
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

