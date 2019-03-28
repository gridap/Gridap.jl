module Polynomials

using Numa
using Numa.Helpers
using Numa.Quadratures
using Base.Cartesian

export MultivariatePolynomialBasis, evaluate!, valuetype, grad

export TensorProductPolynomialBasis
export TensorProductMonomialBasis

export UnivariatePolynomialBasis
export UnivariateMonomialBasis

export derivative, tensorproduct!, gradient, tensorproductsquare!
export gradient

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

valuetype(::Type{R} where R<:MultivariatePolynomialBasis{D,T}) where {D,T} = T

valuetype(self::MultivariatePolynomialBasis) = valuetype(typeof(self))

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

# TODO: This is a temporary dummy implementation that has to be deleted and
# replaced by concrete implementations that use the functionality below.
# It serves now just as an example

export ShapeFunctionsScalarQua4

struct ShapeFunctionsScalarQua4 <: MultivariatePolynomialBasis{2,Float64} end

struct GradShapeFunctionsScalarQua4 <: MultivariatePolynomialBasis{2,VectorValue{2}} end

Base.length(::ShapeFunctionsScalarQua4) = 4

function evaluate!(
  ::ShapeFunctionsScalarQua4,points::Array{Point{2},1},v::Array{Float64,2})
  for (i,point) in enumerate(points)
    xi = point[1]
    eta = point[2]
    v[1,i] = (1-xi)*(1-eta)/4.0
    v[2,i] = (1+xi)*(1-eta)/4.0
    v[3,i] = (1-xi)*(1+eta)/4.0
    v[4,i] = (1+xi)*(1+eta)/4.0
  end
end

gradient(::ShapeFunctionsScalarQua4) = GradShapeFunctionsScalarQua4()

Base.length(::GradShapeFunctionsScalarQua4) = 4

function evaluate!(
  ::GradShapeFunctionsScalarQua4,points::Array{Point{2},1},v::Array{VectorValue{2},2})
  for (i,point) in enumerate(points)
    xi = point[1]
    eta = point[2]
    v[1,i] = VectorValue{2}( (eta-1)/4.0, (xi-1)/4.0 )
    v[2,i] = VectorValue{2}( (1-eta)/4.0,-(1+xi)/4.0 )
    v[3,i] = VectorValue{2}(-(1+eta)/4.0, (1-xi)/4.0 )
    v[4,i] = VectorValue{2}( (1+eta)/4.0, (1+xi)/4.0 )
  end
end

gradient(::GradShapeFunctionsScalarQua4) = @notimplemented

# Previous functionality.
# It has to be used to implement the abstract interface above

end # Polynomials module
