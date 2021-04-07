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
using FillArrays

using Gridap.Helpers
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.ReferenceFEs

import Gridap.ReferenceFEs: num_dims
import Gridap.ReferenceFEs: num_point_dims

export Quadrature
export QuadratureName
export GenericQuadrature
export num_points
export get_coordinates
export get_weights
export get_name
export num_dims
export num_point_dims
export test_quadrature
export tensor_product
export duffy

include("Quadratures.jl")

include("TensorProductQuadratures.jl")

include("DuffyQuadratures.jl")

#include("SymLegendreQuadratures.jl")

end # module
