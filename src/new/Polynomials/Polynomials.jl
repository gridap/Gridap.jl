"""
This module provides a collection of multivariate polynomial bases.

The exported names are:

$(EXPORTS)
"""
module Polynomials

using DocStringExtensions
using LinearAlgebra: mul!
using Gridap.Helpers
using Gridap.Inference
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields

import Gridap.Fields: evaluate_field!
import Gridap.Fields: field_cache
import Gridap.Fields: evaluate_gradient!
import Gridap.Fields: gradient_cache
import Gridap.Fields: evaluate_hessian!
import Gridap.Fields: hessian_cache

export MonomialBasis
export QGradMonomialBasis
export QCurlGradMonomialBasis
export change_basis
export get_exponents

export get_order
export get_orders

include("MonomialBases.jl")

include("QGradMonomialBases.jl")

include("QCurlGradMonomialBases.jl")

include("ChangeBasis.jl")

end # module
