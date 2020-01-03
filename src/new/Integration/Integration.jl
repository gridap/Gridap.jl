"""

The exported names are
$(EXPORTS)
"""
module Integration

using DocStringExtensions
using Test
using Gridap.Helpers
using Gridap.Fields

import Gridap.ReferenceFEs: num_dims
import Gridap.ReferenceFEs: num_point_dims

export Quadrature
export GenericQuadrature
export num_points
export get_coordinates
export get_weights
export test_quadrature

include("Quadratures.jl")

end # module
