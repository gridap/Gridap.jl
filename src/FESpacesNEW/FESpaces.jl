"""

The exported names are
$(EXPORTS)
"""
module FESpaces

using DocStringExtensions
using Test
using FillArrays
using SparseArrays

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData

import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Geometry: get_triangulation
import Gridap.Geometry: get_cell_shapefuns
import Gridap.CellData: attach_constraints_rows
import Gridap.CellData: attach_constraints_cols
import Gridap.CellData: CellField
import Gridap.CellData: get_cell_data
import Gridap.CellData: DomainStyle

export FEFunction
export get_free_values
export get_cell_dof_values
export get_fe_space
export test_fe_function

export FESpace
export ConstraintStyle
export num_free_dofs
export zero_free_values
export EvaluationFunction
export get_cell_dof_ids
export get_cell_shapefuns
export get_cell_shapefuns_trial
export has_constraints
export get_cell_constraints
export get_cell_isconstrained
export get_dof_value_type
export test_fe_space

export SingleFieldFESpace
export SingleFieldFEFunction
export get_cell_dof_basis
export num_dirichlet_dofs
export zero_dirichlet_values
export num_dirichlet_tags
export get_dirichlet_dof_tag
export scatter_free_and_dirichlet_values
export gather_free_and_dirichlet_values!
export gather_dirichlet_values
export gather_dirichlet_values!
export gather_free_values
export gather_free_values!
export test_single_field_fe_space
export get_dirichlet_values
export interpolate
export interpolate!
export interpolate_everywhere
export interpolate_everywhere!
export interpolate_dirichlet
export interpolate_dirichlet!
export compute_dirichlet_values_for_tags
export compute_dirichlet_values_for_tags!

export UnconstrainedFESpace
export compute_conforming_cell_dofs
export compute_cell_space

export TrialFESpace
export TrialFESpace!
export HomogeneousTrialFESpace
export HomogeneousTrialFESpace!

include("FESpaceInterface.jl")

include("SingleFieldFESpaces.jl")

include("UnconstrainedFESpaces.jl")

include("ConformingFESpaces.jl")

include("FESpaceFactories.jl")

include("TrialFESpaces.jl")

end # module
