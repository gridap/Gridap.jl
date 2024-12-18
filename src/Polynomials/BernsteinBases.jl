"""
    Bernstein <: Polynomial

Type representing Bernstein polynomials
"""
struct Bernstein <: Polynomial end

isHierarchical(::Type{Bernstein}) = false

"""
    BernsteinBasis{D,V,K} = TensorPolynomialBasis{D,V,K,Bernstein}

Multivariate scalar' or `Multivalue`'d monomial basis, see [`TensorPolynomialBasis`](@ref)
"""
const BernsteinBasis{D,V,K} = TensorPolynomialBasis{D,V,K,Bernstein}

"""
    BernsteinBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector)
    BernsteinBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function])
    BernsteinBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function])

Convenience constructors of [`BernsteinBasis`](@ref).
"""
BernsteinBasis(args...) = TensorPolynomialBasis(Bernstein, args...)

QGradBernsteinBasis(args...)     = QGradBasis(Bernstein, args...)
QCurlGradBernsteinBasis(args...) = QCurlGradBasis(Bernstein, args...)


# 1D evaluation implementation

# TODO Optimize with in-place De Casteljau

"""
    binoms(::Val{K})

Returns the tuple of binomials ( C₍ₖ₀₎, C₍ₖ₁₎, ..., C₍ₖₖ₎ )
"""
binoms(::Val{K}) where K = ntuple( i -> binomial(K,i-1), Val(K+1))


# jth Bernstein poly of order K at x:
# Bᵏⱼ(x) = binom(K,j) * x^j * (1-x)^(K-j)
function _evaluate_1d!(::Type{Bernstein},::Val{K}, v::AbstractMatrix{T},x,d) where {K,T<:Number}
  b = binoms(Val(K))
  n = K + 1
  @inbounds λ1 = x[d]
  λ2 = one(T) - λ1

  for i in 1:n
    j = i-1
    λ1_j = λ1^(j)
    λ2_j = λ2^(K-j)
    @inbounds v[d,i] = b[i] * λ1_j * λ2_j # order K
  end
end


function _gradient_1d!(::Type{Bernstein},::Val{0},g::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds g[d,1] = zero(T)
end
function _gradient_1d!(::Type{Bernstein},::Val{1},g::AbstractMatrix{T},x,d) where {T<:Number}
  o = one(T)
  @inbounds g[d,1] = -o
  @inbounds g[d,2] =  o
end

# First derivative of the jth Bernstein poly of order K at x:
# (Bᵏⱼ)'(x) = K * ( Bᵏ⁻¹ⱼ₋₁(x) - Bᵏ⁻¹ⱼ(x) )
#           = K * x^(j-1) * (1-x)^(K-j-1) * ((1-x)*binom(K-1,j-1) - x*binom(K-1,j))
function _gradient_1d!(::Type{Bernstein},::Val{K}, g::AbstractMatrix{T},x,d) where {K,T<:Number}
  b = binoms(Val(K-1))
  n = K + 1

  @inbounds λ1 = x[d]
  λ2 = one(T) - λ1

  @inbounds g[1] =-K * λ2^(K-1)
  @inbounds g[n] = K * λ1^(K-1)
  for i in 2:n-1
    j = i-1
    λ1_j = λ1^(j-1)
    λ2_j = λ2^(K-j-1)
    @inbounds g[d,i] = K * λ1_j * λ2_j *(λ2*b[i-1] - λ1*b[i]) # order K-1
  end
end


function _hessian_1d!(::Type{Bernstein},::Val{0},h::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds h[d,1] = zero(T)
end
function _hessian_1d!(::Type{Bernstein},::Val{1},h::AbstractMatrix{T},x,d) where {T<:Number}
  @inbounds h[d,1] = zero(T)
  @inbounds h[d,2] = zero(T)
end
function _hessian_1d!(::Type{Bernstein},::Val{2},h::AbstractMatrix{T},x,d) where {T<:Number}
  o = one(T)
  @inbounds h[d,1] =  2o
  @inbounds h[d,2] = -4o
  @inbounds h[d,3] =  2o
end

# Second derivative of the jth Bernstein poly of order K at x:
# (Bᵏⱼ)''(x) = K(K-1) * ( Bᵏ⁻²ⱼ₋₂(x) -2*Bᵏ⁻²ⱼ₋₁(x) + Bᵏ⁻²ⱼ(x) )
#            = K(K-1) * x^(j-2) * (1-x)^(K-j-2) * ( (1-x)^2*binom(K-2,j-2)
#                  - 2x*(1-x)*binom(K-2,j-1) + (x)^2*binom(K-2,j)
#              )
function _hessian_1d!(::Type{Bernstein},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number}
  b = binoms(Val(K-2))
  n = K + 1
  C = K*(K-1)

  @inbounds λ1 = x[d]
  λ2 = one(T) - λ1

  @inbounds h[1] = C *    λ2^(K-2)
  @inbounds h[2] = C * (-2λ2^(K-2) + (K-2)*λ2^(K-3)*λ1)
  @inbounds h[n-1]=C * (-2λ1^(K-2) + (K-2)*λ1^(K-3)*λ2)
  @inbounds h[n] = C *    λ1^(K-2)
  for i in 3:n-2
    j = i-1
    λ1_j = λ1^(j-2)
    λ2_j = λ2^(K-j-2)
    @inbounds h[d,i] = C * λ1_j * λ2_j *(λ2*λ2*b[i-2] -2λ2*λ1*b[i-1] + λ1*λ1*b[i]) # order K-2
  end
end
