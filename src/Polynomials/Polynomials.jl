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
export isHierarchical
export PolynomialBasis

export TensorPolynomialBasis
export get_exponents
export get_order
export get_orders
export num_terms

export MonomialBasis
export LegendreBasis
export ChebyshevBasis
export ModalC0Basis
export BernsteinBasis

export QGradMonomialBasis
export QGradLegendreBasis
export QGradChebyshevBasis

export QCurlGradMonomialBasis
export QCurlGradLegendreBasis
export QCurlGradChebyshevBasis

export PCurlGradMonomialBasis
export PCurlGradLegendreBasis


include("PolynomialInterfaces.jl")

include("TensorPolynomialBases.jl")

include("CompWiseTensorPolyBases.jl")

include("NonTensorRTPolyBasis.jl")

include("MonomialBases.jl")

include("LegendreBases.jl")

#include("JacobiBases.jl")

include("ChebyshevBases.jl")

include("BernsteinBases.jl")

include("ModalC0Bases.jl")

include("NedelecPrebasisOnSimplex.jl")

end # module
