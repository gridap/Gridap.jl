#################################################
# Polynomial family/types with their parameters #
#################################################

"""
Abstract type for 1- or n-D polynomial family of maximum order K.

`PolynomialType` subtypes must keep `K` as last type parameter.

Implement [`isHierarchical`](@ref) and [`get_order`](@ref).
"""
abstract type PolynomialType{K} end

"""
Return true if the polynomial basis of order `K` of the given type is the union
of the basis of order `K-1` and an other order `K` polynomial.
"""
isHierarchical(::PolynomialType) = @abstractmethod

"""
    get_order(PT<:PolynomialType{K}}) returns K
"""
get_order(::Type{PolynomialType{K}}) where K = K


#"""
#    JacobiPType{α, β, K} <: PolynomialType{K}
#
#Type representing Jacobi polynomials of parameters `α` and `β` and order up to `K`
#"""
#struct JacobiPType{α, β, K} <: PolynomialType{K} end
#isHierarchical(::JacobiPType) = true


######################
# 1D polynomial APIs #
######################

"""
    _evaluate_1d!(::Type{<:PolynomialType{K}}, v, x, d)

Evaluates in place the 1D basis polynomials of the given type at one nD point `x`
along the given coordinate 1 ≤ `d` ≤ nD.

`v` is an AbstractMatrix of size (at least} d×(K+1), such that the 1 ≤ i ≤ `K`+1
values are stored in `v[d,i]`.
"""
function _evaluate_1d!(
  ::Type{<:PolynomialType{K}}, v::AbstractMatrix{T},x,d) where {K,T<:Number}

  @abstractmethod
end

"""
    _gradient_1d!(::Type{<:PolynomialType{K}}, g, x, d)

Like [`_evaluate_1d!`](@ref), but computes the first derivative of the basis functions.
"""
function _gradient_1d!(
  ::Type{<:PolynomialType{K}}, g::AbstractMatrix{T},x,d) where {K,T<:Number}

  @abstractmethod
end

"""
    _hessian_1d!(::Type{<:PolynomialType{K}}, h, x, d)

Like [`_evaluate_1d!`](@ref), but computes the second derivative of the basis functions.
"""
function _hessian_1d!(
  ::Type{<:PolynomialType{K}}, h::AbstractMatrix{T},x,d) where {K,T<:Number}

  @abstractmethod
end

"""
    _derivatives_1d!(PT::Type{<:PolynomialType{K}}, (v,g,...), x, d)

Same as calling
```
_evaluate_1d!(PT, v, x d)
_gradient_1d!(PT, g, x d)
          ⋮
```
but with possible performence optimization.
"""
function _derivatives_1d!(PT::Type{<:PolynomialType},t::NTuple{N},x,d) where N
  @abstractmethod
end

function _derivatives_1d!(PT::Type{<:PolynomialType},t::NTuple{1},x,d)
  _evaluate_1d!(PT, t[1], x, d)
end

function _derivatives_1d!(PT::Type{<:PolynomialType},t::NTuple{2},x,d)
  _evaluate_1d!(PT, t[1], x, d)
  _gradient_1d!(PT, t[2], x, d)
end

function _derivatives_1d!(PT::Type{<:PolynomialType},t::NTuple{3},x,d)
  _evaluate_1d!(PT, t[1], x, d)
  _gradient_1d!(PT, t[2], x, d)
  _hessian_1d!( PT, t[3], x, d)
end
