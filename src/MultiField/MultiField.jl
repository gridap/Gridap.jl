"""

The exported names are
$(EXPORTS)
"""
module MultiField

using DocStringExtensions
using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.Integration
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields

using Gridap.FESpaces: FEBasis, TestBasis, TrialBasis, get_cell_dof_values
using Gridap.Arrays: BlockArrayCooMap
using Gridap.ReferenceFEs: HEX_AXIS, TET_AXIS,
                           get_faces, get_grid_topology, get_node_coordinates,
                           num_cell_dims

using FillArrays
using SparseArrays
using LinearAlgebra
using BlockArrays
using NearestNeighbors
using StaticArrays

import Gridap.Arrays: evaluate!
import Gridap.Arrays: return_cache

export num_fields
export compute_field_offsets
export restrict_to_field
export MultiFieldCellField
export MultiFieldFESpace
export MultiFieldFEFunction
export MultiFieldStyle
export ConsecutiveMultiFieldStyle

include("MultiFieldCellFields.jl")

include("MultiFieldFESpaces.jl")

include("MultiFieldFEFunctions.jl")

end # module
