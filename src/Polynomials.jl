module Polynomials

using Numa #@fverdugo to be eliminated
using Numa.Helpers
using Numa.FieldValues
using Numa.Quadratures

using Base.Cartesian

export MultivariatePolynomialBasis
export TensorProductPolynomialBasis
export TensorProductMonomialBasis
export UnivariatePolynomialBasis
export UnivariateMonomialBasis

export evaluate!
export gradient

# @fverdugo: really needed to export?
# If they are just needed in the tests use qualified names there
export derivative, tensorproduct!, tensorproductsquare!

# @santiagobadia : To put my polynomials in accordance with Abstract methods
# Abstract types and interfaces

"""
Abstract type representing a multivariate polynomial basis
with value of type T in a coordinate space of D dimensions
"""
abstract type MultivariatePolynomialBasis{D,T} end

Base.length(::MultivariatePolynomialBasis)::Int = @abstractmethod

"""
First axis of v for dofs, second for points
"""
evaluate!(::MultivariatePolynomialBasis{D,T},::Array{Point{D},1},v::Array{T,2}) where {D,T} = @abstractmethod

"""
Returns a MultivariatePolynomialBasis{TG,D} where TG
is a type whose rank is one unit grater than the one of T
"""
gradient(::MultivariatePolynomialBasis{D,T} where{D,T})::MultivariatePolynomialBasis{D,TG} = @abstractmethod

"""
Abstract basis of univariate polynomials
"""
abstract type UnivariatePolynomialBasis end

# Concrete structs

"""
Univariate monomial basis of a given `order`
"""
struct UnivariateMonomialBasis <: UnivariatePolynomialBasis
  order::Int64
end

"""
Multivariate polynomial basis obtained as tensor product of univariate polynomial basis
per dimension
"""
struct TensorProductPolynomialBasis
  polynomials::Vector{UnivariatePolynomialBasis}
end

# Methods

include("PolynomialsMethods.jl")

end # Polynomials module
