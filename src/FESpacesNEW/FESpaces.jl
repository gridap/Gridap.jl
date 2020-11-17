"""

The exported names are
$(EXPORTS)
"""
module FESpaces

using DocStringExtensions
using Test
using FillArrays
using BlockArrays
using SparseArrays
using LinearAlgebra

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.TensorValues

import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!
import Gridap.Geometry: get_triangulation
import Gridap.Geometry: get_cell_shapefuns
import Gridap.CellData: attach_constraints_rows
import Gridap.CellData: attach_constraints_cols
import Gridap.CellData: CellField
import Gridap.CellData: get_cell_data
import Gridap.CellData: DomainStyle
import Gridap.CellData: change_domain_skeleton
import Gridap.CellData: gradient
import Gridap.CellData: ∇∇

import Gridap.Algebra: allocate_residual
import Gridap.Algebra: allocate_jacobian
import Gridap.Algebra: residual!
import Gridap.Algebra: jacobian!
import Gridap.Algebra: residual
import Gridap.Algebra: jacobian
import Gridap.Algebra: residual_and_jacobian!
import Gridap.Algebra: residual_and_jacobian
import Gridap.Algebra: zero_initial_guess
import Gridap.Algebra: get_matrix
import Gridap.Algebra: get_vector
import Gridap.Algebra: solve!
import Gridap.Algebra: solve
import Gridap.Algebra: allocate_vector
import Gridap.Algebra: allocate_matrix
import Gridap.Algebra: allocate_matrix_and_vector

export FEFunction
export get_free_values
export get_cell_dof_values
export get_fe_space
export test_fe_function

export FESpace
export ConstraintStyle
export Constrained
export UnConstrained
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

export TestFESpace
export TrialFESpace
export TrialFESpace!
export HomogeneousTrialFESpace
export HomogeneousTrialFESpace!

export Assembler
export AssemblyStrategy
export row_map
export col_map
export row_mask
export col_mask
export DefaultAssemblyStrategy
export get_test
export get_trial
export assemble_matrix!
export assemble_matrix_add!
export assemble_matrix
export assemble_vector!
export assemble_vector_add!
export assemble_vector
export assemble_matrix_and_vector!
export assemble_matrix_and_vector_add!
export assemble_matrix_and_vector
export allocate_vector
export allocate_matrix
export allocate_matrix_and_vector
export test_assembler
export collect_cell_matrix
export collect_cell_vector
export collect_cell_matrix_and_vector
export get_matrix_type
export get_vector_type
export SparseMatrixAssembler
export GenericSparseMatrixAssembler
export count_matrix_nnz_coo
export count_matrix_and_vector_nnz_coo
export fill_matrix_coo_symbolic!
export fill_matrix_and_vector_coo_symbolic!
export fill_matrix_coo_numeric!
export fill_matrix_and_vector_coo_numeric!
export test_sparse_matrix_assembler

export FEOperator
export AffineFEOperator
export test_fe_operator
export get_algebraic_operator

export FESolver
export test_fe_solver
export LinearFESolver
export NonlinearFESolver

export FESpaceWithConstantFixed
export ZeroMeanFESpace
export CLagrangianFESpace
export DirichletFESpace
export FESpaceWithLinearConstraints

include("FESpaceInterface.jl")

include("SingleFieldFESpaces.jl")

include("UnconstrainedFESpaces.jl")

include("ConformingFESpaces.jl")

include("FESpaceFactories.jl")

include("TrialFESpaces.jl")

include("Assemblers.jl")

include("SparseMatrixAssemblers.jl")

include("FEOperators.jl")

include("AffineFEOperators.jl")

include("FEOperatorsFromWeakForm.jl")

include("FESolvers.jl")

include("DiscontinuousFESpaces.jl")

include("FESpacesWithConstantFixed.jl")

include("ZeroMeanFESpaces.jl")

include("CLagrangianFESpaces.jl")

include("DirichletFESpaces.jl")

include("ExtendedFESpaces.jl")

include("FESpacesWithLinearConstraints.jl")

end # module
