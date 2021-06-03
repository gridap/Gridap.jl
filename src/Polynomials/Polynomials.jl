"""
This module provides a collection of multivariate polynomial bases.

The exported names are:

$(EXPORTS)
"""
module Polynomials

using DocStringExtensions
using LinearAlgebra: mul!
using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields

import Gridap.Fields: evaluate!
import Gridap.Fields: return_cache
import Gridap.Arrays: return_type

export MonomialBasis
export QGradMonomialBasis
export QCurlGradMonomialBasis
export PCurlGradMonomialBasis
export get_exponents

export get_order
export get_orders
export num_terms

include("MonomialBases.jl")

include("QGradMonomialBases.jl")

include("QCurlGradMonomialBases.jl")

include("PCurlGradMonomialBases.jl")

end # module
