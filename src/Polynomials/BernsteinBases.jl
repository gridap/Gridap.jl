"""
    Bernstein <: Polynomial

Type representing Bernstein polynomials
"""
struct Bernstein <: Polynomial end

isHierarchical(::Bernstein) = false

"""
    BernsteinBasis{D,V,K} = TensorPolynomialBasis{D,V,K,Bernstein}

Multivariate scalar' or `Multivalue`'d monomial basis, see [`TensorPolynomialBasis`](@ref)
"""
const BernsteinBasis{D,V,K} = TensorPolynomialBasis{D,V,K,Bernstein}

"""
    BernsteinBasis{D}(::Type{V}, order::Int, terms::Vector) where {D,V}
    BernsteinBasis{D}(::Type{V}, orders::Tuple [, filter::Function]) where {D,V}
    BernsteinBasis{D}(::Type{V}, order::Int [, filter::Function]) where {D,V}

Convenience constructors of BernsteinBasis{D,V}.
"""
BernsteinBasis{D}(args...) where {D} = TensorPolynomialBasis{D}(Bernstein, args...)

# 1D evaluation implementation

# TODO Optimize with in-place De Casteljau

binoms(::Val{K}) where K = ntuple( i -> binomial(K,i-1), K+1)

# jth Bernstein poly of order k at x:
# Bᵏⱼ(x) = binom(k,j) * x^j * (1-x)^(k-j)
function _evaluate_1d!(::Type{Bernstein}, k, v::AbstractMatrix{T},x,d) where T<:Number
  b = binoms(Val(k))
  n = k + 1
  @inbounds λ1 = x[d]
  λ2 = one(T) - λ1
  for i in 1:n
    j = i-1
    λ1_j = λ1^(j)
    λ2_j = λ2^(k-j)
    @inbounds v[d,i] = b[i] * λ1_j * λ2_j # order k
  end
end

# First derivative of the jth Bernstein poly of order k at x:
# (Bᵏⱼ)'(x) = k * ( Bᵏ⁻¹ⱼ₋₁(x) - Bᵏ⁻¹ⱼ(x) )
#           = k * x^(j-1) * (1-x)^(k-j-1) * ((1-x)*binom(k-1,j-1) - x*binom(k-1,j))
function _gradient_1d!(::Type{Bernstein}, k, g::AbstractMatrix{T},x,d) where T<:Number
  if iszero(k)
    @inbounds g[d,1] = zero(T)
    return
  end

  o = one(T)
  if isone(k)
    @inbounds g[d,1] = -o
    @inbounds g[d,2] =  o
    return
  end

  b = binoms(Val(k-1))
  n = k + 1
  @inbounds λ1 = x[d]
  λ2 = o - λ1

  @inbounds g[1] =-k * λ2^(k-1)
  @inbounds g[n] = k * λ1^(k-1)
  for i in 2:n-1
    j = i-1
    λ1_j = λ1^(j-1)
    λ2_j = λ2^(k-j-1)
    @inbounds g[d,i] = k * λ1_j * λ2_j *(λ2*b[i-1] - λ1*b[i]) # order k-1
  end
end

# Second derivative of the jth Bernstein poly of order k at x:
# (Bᵏⱼ)''(x) = k(k-1) * ( Bᵏ⁻²ⱼ₋₂(x) -2*Bᵏ⁻²ⱼ₋₁(x) + Bᵏ⁻²ⱼ(x) )
#            = k(k-1) * x^(j-2) * (1-x)^(k-j-2) * ( (1-x)^2*binom(k-2,j-2)
#                  - 2x*(1-x)*binom(k-2,j-1) + (x)^2*binom(k-2,j)
#              )
function _hessian_1d!(::Type{Bernstein}, k, h::AbstractMatrix{T},x,d) where T<:Number
  z = zero(T)
  if iszero(k)
    @inbounds h[d,1] = z
    return
  end
  if isone(k)
    @inbounds h[d,1] = z
    @inbounds h[d,2] = z
    return
  end

  o = one(T)
  if k == 2
    @inbounds h[d,1] =  2o
    @inbounds h[d,2] = -4o
    @inbounds h[d,3] =  2o
    return
  end

  b = binoms(Val(k-2))
  n = k + 1
  C = k*(k-1)
  @inbounds λ1 = x[d]
  λ2 = o - λ1

  @inbounds h[1] = C *    λ2^(k-2)
  @inbounds h[2] = C * (-2λ2^(k-2) + (k-2)*λ2^(k-3)*λ1)
  @inbounds h[n-1]=C * (-2λ1^(k-2) + (k-2)*λ1^(k-3)*λ2)
  @inbounds h[n] = C *    λ1^(k-2)
  for i in 3:n-2
    j = i-1
    λ1_j = λ1^(j-2)
    λ2_j = λ2^(k-j-2)
    @inbounds h[d,i] = C * λ1_j * λ2_j *(λ2*λ2*b[i-2] -2λ2*λ1*b[i-1] + λ1*λ1*b[i]) # order k-2
  end
end
