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
export Monomial
export Legendre
export Chebyshev
export ModalC0
export Bernstein

export PolynomialBasis
export get_order

export UniformPolyBasis
export get_exponents
export get_orders
export MonomialBasis
export LegendreBasis
export ChebyshevBasis
export BernsteinBasis

export CompWiseTensorPolyBasis
export QGradBasis
export QCurlGradBasis

export NedelecPolyBasisOnSimplex
export PGradBasis

export RaviartThomasPolyBasis
export PCurlGradBasis

export ModalC0Basis

# deprecated
export num_terms
export QGradMonomialBasis
export QCurlGradMonomialBasis
export PCurlGradMonomialBasis


include("PolynomialInterfaces.jl")

include("UniformPolyBases.jl")

include("CompWiseTensorPolyBases.jl")

include("NedelecPolyBases.jl")

include("RaviartThomasPolyBases.jl")

include("MonomialBases.jl")

include("LegendreBases.jl")

include("ChebyshevBases.jl")

include("BernsteinBases.jl")

include("ModalC0Bases.jl")

include("Deprecated.jl")

end # module
