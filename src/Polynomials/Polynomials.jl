"""
This module provides a collection of multivariate polynomial bases.

The exported names are:

$(EXPORTS)
"""
module Polynomials

using DocStringExtensions
using LinearAlgebra: mul!
using StaticArrays
using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields

using PolynomialBases: jacobi, jacobi_and_derivative

import Gridap.Fields: evaluate!
import Gridap.Fields: return_cache
import Gridap.Arrays: return_type

export Polynomial
export MonomialBasis
export QGradMonomialBasis
export QCurlGradMonomialBasis
export PCurlGradMonomialBasis
export ModalC0Basis
export LegendreBasis
export QGradLegendrePolynomialBasis
export QCurlGradLegendrePolynomialBasis
export PCurlGradLegendrePolynomialBasis
export ChebyshevBasis
export QGradChebyshevPolynomialBasis
export QCurlGradChebyshevPolynomialBasis
export BernsteinBasis

export get_exponents
export get_order
export get_orders
export num_terms
export isHierarchical

include("PolynomialInterfaces.jl")

include("TensorPolynomialBases.jl")

include("MonomialBases.jl")

include("OldMonomialHelpers.jl")

include("QGradMonomialBases.jl")

include("QCurlGradMonomialBases.jl")

include("PCurlGradMonomialBases.jl")

include("ModalC0Bases.jl")

include("LegendreBases.jl")

include("QGradLegendrePolynomialBases.jl")

include("PCurlGradLegendrePolynomialBases.jl")

include("ChebyshevBases.jl")

include("BernsteinBases.jl")

end # module
