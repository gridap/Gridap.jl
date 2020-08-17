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
using Gridap.Polynomials
using Gridap.Integration
using Gridap.ReferenceFEs

import Gridap.Arrays: get_array
import Gridap.Fields: evaluate
import Gridap.Fields: gradient
import Gridap.Helpers: operate

export CellField
export get_cell_map
export get_cell_axes
export RefStyle
export is_in_ref_space
export is_in_physical_space
export MetaSizeStyle
export get_metasize
export is_basis
export is_test
export is_trial
export test_cell_field
export GenericCellField
export similar_object
export change_ref_style
export to_ref_space
export to_physical_space
export trialize_cell_basis
export convert_to_cell_field
export SkeletonCellField
export get_inward
export get_outward
export jump
export mean
export merge_cell_dof_ids_at_skeleton
export merge_cell_fields_at_skeleton

include("CellFields.jl")

end # module

