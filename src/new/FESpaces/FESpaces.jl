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
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.Integration
using Gridap.Algebra

using Gridap.Geometry: UnimplementedField

import Gridap.Arrays: get_array
import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Geometry: get_cell_map
import Gridap.Geometry: get_cell_shapefuns
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_cell_type
import Gridap.Geometry: operate_cell_field
import Gridap.Geometry: similar_cell_field
import Gridap.Geometry: jump
import Gridap.Geometry: mean
import Gridap.Geometry: restrict
import Gridap.Fields: integrate

import Gridap.Algebra: allocate_residual
import Gridap.Algebra: allocate_jacobian
import Gridap.Algebra: residual!
import Gridap.Algebra: jacobian!
import Gridap.Algebra: residual
import Gridap.Algebra: jacobian
import Gridap.Algebra: zero_initial_guess
import Gridap.Algebra: get_matrix
import Gridap.Algebra: get_vector

export FEFunctionStyle
export is_a_fe_function
export get_free_values
export get_fe_space
export test_fe_function

export FESpace
export FEFunction
export num_free_dofs
export get_cell_basis
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
export num_dirichlet_tags
export gather_free_values
export get_dirichlet_dof_tag
export compute_free_and_dirichlet_values
export compute_dirichlet_values
export compute_free_values
export inpterpolate
export compute_dirichlet_values_for_tags
export test_single_field_fe_space
export interpolate
export interpolate_everywhere
export interpolate_dirichlet
export get_cell_dof_basis

export SingleFieldFEFunction

export UnsconstrainedFESpace
export GradConformingFESpace

export CellBasis
export test_cell_basis
export CellMatrixField
export test_cell_matrix_field
export GenericCellBasis
export GenericCellMatrixField
export TrialStyle
export is_trial
export is_test

export TrialFESpace
export compute_conforming_cell_dofs
export SparseMatrixAssembler

export FEOperator
export test_fe_operator
export AffineFEOperator

include("CellBases.jl")

include("FEFunctions.jl")

include("FESpacesInterfaces.jl")

include("Assemblers.jl")

include("FEOperators.jl")

include("SingleFieldFESpaces.jl")

include("SingleFieldFEFunctions.jl")

include("TrialFESpaces.jl")

include("SparseMatrixAssemblers.jl")

include("UnconstrainedFESpaces.jl")

include("ConformingFESpaces.jl")

include("AffineFEOperators.jl")

end # module
