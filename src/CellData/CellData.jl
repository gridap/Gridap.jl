"""

The exported names are
$(EXPORTS)
"""
module CellData

using Test
using DocStringExtensions
using FillArrays

using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Geometry
using Gridap.Integration

import Gridap.Arrays: get_array
import Gridap.Arrays: evaluate!
import Gridap.Integration: get_coordinates
import Gridap.Integration: get_weights

export DomainStyle
export ReferenceDomain
export PhysicalDomain
export CellDatum
export get_cell_data
export get_triangulation
export change_domain
export test_cell_datum
export CellPoint
export get_cell_points
export CellField
export GenericCellField
export CellWeight
export CellQuadrature

include("CellDataInterface.jl")

include("CellFields.jl")

include("CellQuadratures.jl")

end # module
