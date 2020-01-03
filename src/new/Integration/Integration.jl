"""

The exported names are
$(EXPORTS)
"""
module Integration

using DocStringExtensions
using Test
using QuadGK: gauss
using FastGaussQuadrature: gaussjacobi
using FastGaussQuadrature: gausslegendre
using Gridap.Helpers
using Gridap.Fields
using Gridap.TensorValues

import Gridap.ReferenceFEs: num_dims
import Gridap.ReferenceFEs: num_point_dims

export Quadrature
export GenericQuadrature
export TensorProductQuadrature
export DuffyQuadrature
export num_points
export get_coordinates
export get_weights
export test_quadrature

include("Quadratures.jl")

include("TensorProductQuadratures.jl")

include("DuffyQuadratures.jl")

end # module
