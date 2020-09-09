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
using Gridap.Algebra

using Gridap.Arrays: IdentityVector

import Gridap.Arrays: get_array
import Gridap.Arrays: lazy_append
import Gridap.Fields: evaluate
import Gridap.Fields: gradient
import Gridap.Fields: integrate
import Gridap.Helpers: operate
import Gridap.Integration: get_coordinates
import Gridap.Integration: get_weights
import Gridap.Fields: field_cache
import Gridap.Fields: evaluate_field!
import Gridap.Fields: evaluate_field_array
import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Arrays: reindex

export CellPoint
export GenericCellPoint
export get_cell_id
export num_cell_ids
export test_cell_point

export CellMap
export GenericCellMap
export test_cell_map

export CellField
export get_cell_axes
export get_memo
export MetaSizeStyle
export get_metasize
export is_basis
export is_test
export is_trial
export test_cell_field
export GenericCellField
export trialize_cell_basis
export convert_to_cell_field
export get_inward
export get_outward
export jump
export mean
export merge_cell_dofs_at_skeleton
export merge_cell_fields_at_skeleton
export inverse_map
export SkeletonFaceMap
export FaceMap
export ReindexedCellMap

export CellQuadrature
export get_coordinates
export get_weights

export QPointCellField
export update_state_variables!

export CellDof
export test_cell_dof_basis
export GenericCellDof
export MappedCellDof

export attach_dirichlet

export attach_constraints_rows
export attach_constraints_cols
export merge_cell_constraints_at_skeleton
export identity_constraints

export @law
export operate
export GridapType
export ∫

include("CellPoints.jl")

include("CellMaps.jl")

include("CellFields.jl")

include("CellQuadratures.jl")

#include("QPointCellFields.jl")
#
#include("CellDofs.jl")
#
#include("AttachDirichlet.jl")
#
#include("AttachConstraints.jl")
#
#include("Law.jl")

end # module

