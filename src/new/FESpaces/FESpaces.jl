"""

The exported names are
$(EXPORTS)
"""
module FESpaces

using DocStringExtensions
using Test

using Gridap.Helpers
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Fields

import Gridap.Arrays: get_array
import Gridap.Geometry: get_cell_map
import Gridap.Geometry: get_cell_shapefuns
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_cell_type

export FEFunction
export get_free_values
export get_fe_space
export test_fe_function

export FESpace
export num_free_dofs
export get_cell_fe_basis
export zero_free_values
export apply_constraints_matrix_cols
export apply_constraints_matrix_rows
export apply_constraints_vector
export apply_constraints_matrix_and_vector_rows
export test_fe_space

export Assembler
export get_test
export get_trial
export allocate_matrix
export assemble_matrix!
export assemble_matrix
export allocate_vector
export assemble_vector!
export assemble_vector
export allocate_matrix_and_vector
export assemble_matrix_and_vector!
export assemble_matrix_and_vector
export test_assembler

export SingleFieldFESpace
export num_dirichlet_dofs
export get_cell_fe
export get_cell_dofs
export zero_dirichlet_values
export gather_free_and_dirichlet_values
export scatter_free_and_dirichlet_values
export get_dirichlet_values
export gather_dirichlet_values
export gather_free_values
export compute_free_and_dirichlet_values
export compute_dirichlet_values
export compute_free_values
export inpterpolate
export compute_dirichlet_values_for_tags
export test_single_field_fe_space

export SingleFieldFEFunction
export SingleFieldCellBasis

export SingleFieldCellFE
export get_cell_dof_basis
export test_single_field_cell_fe

export RefFEsWithMap

include("FEFunctions.jl")

include("FESpacesInterfaces.jl")

include("Assemblers.jl")

include("SingleFieldFESpaces.jl")

end # module
