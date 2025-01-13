"""

$(public_names_in_md(@__MODULE__))
"""
module FESpaces

using DocStringExtensions
using Test
using FillArrays
using SparseArrays
using LinearAlgebra
using StaticArrays
using ForwardDiff

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.TensorValues

using Gridap.Arrays: Reindex, ConfigMap, DualizeMap, AutoDiffMap, lazy_map

using Gridap.Fields: ArrayBlock, BlockMap

using Gridap.Geometry: SkeletonPair

using Gridap.CellData: SkeletonCellFieldPair

import Gridap.Fields: gradient
import Gridap.Fields: ∇∇
import Gridap.Fields: DIV
import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!
import Gridap.Arrays: autodiff_array_gradient # overloaded for Skeleton terms
import Gridap.Arrays: autodiff_array_jacobian
import Gridap.Arrays: autodiff_array_hessian
import Gridap.Geometry: get_triangulation
import Gridap.Geometry: get_cell_shapefuns
import Gridap.Geometry: get_cell_type
import Gridap.Geometry: MappedGrid
import Gridap.Geometry: MappedDiscreteModel
import Gridap.Geometry: get_node_coordinates
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_cell_node_ids
import Gridap.Geometry: OrientationStyle
import Gridap.Geometry: RegularityStyle
import Gridap.Geometry: get_facet_normal
import Gridap.CellData: attach_constraints_rows
import Gridap.CellData: attach_constraints_cols
import Gridap.CellData: CellField
import Gridap.CellData: get_data
import Gridap.CellData: DomainStyle
import Gridap.CellData: change_domain
import Gridap.CellData: change_domain_ref_ref
import Gridap.CellData: change_domain_phys_phys

import Gridap.Algebra: allocate_residual
import Gridap.Algebra: allocate_jacobian
import Gridap.Algebra: residual!
import Gridap.Algebra: jacobian!
import Gridap.Algebra: residual
import Gridap.Algebra: jacobian
import Gridap.Algebra: hessian
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
import Gridap.ReferenceFEs: get_shapefuns
import Gridap.ReferenceFEs: get_dof_basis

export FEFunction
export get_free_dof_values
export get_cell_dof_values
export get_fe_space
export test_fe_function

export FESpace
export ConstraintStyle
export Constrained
export UnConstrained
export num_free_dofs
export get_free_dof_ids
export zero_free_values
export EvaluationFunction
export get_cell_dof_ids
export get_fe_basis
export get_trial_fe_basis
export has_constraints
export get_cell_constraints
export get_cell_isconstrained
export get_cell_is_dirichlet
export get_dof_value_type
export test_fe_space

export SingleFieldFESpace
export SingleFieldFEFunction
export get_fe_dof_basis
export num_dirichlet_dofs
export get_dirichlet_dof_ids
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
export get_dirichlet_dof_values
export interpolate
export interpolate!
export interpolate_everywhere
export interpolate_everywhere!
export interpolate_dirichlet
export interpolate_dirichlet!
export compute_dirichlet_values_for_tags
export compute_dirichlet_values_for_tags!

export UnconstrainedFESpace
export CellConformity
export CellFE
export compute_conforming_cell_dofs
export compute_cell_space

export TestFESpace
export TrialFESpace
export TrialFESpace!
export HomogeneousTrialFESpace
export HomogeneousTrialFESpace!

export FEBasis
export BasisStyle

export Assembler
export AssemblyStrategy
export row_map
export col_map
export row_mask
export col_mask
export DefaultAssemblyStrategy
export GenericAssemblyStrategy
export get_test
export get_trial
export get_rows
export get_cols
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
export get_matrix_builder
export get_vector_builder
export get_matrix_type
export get_vector_type
export SparseMatrixAssembler
export GenericSparseMatrixAssembler
export symbolic_loop_matrix!
export symbolic_loop_vector!
export symbolic_loop_matrix_and_vector!
export numeric_loop_matrix!
export numeric_loop_vector!
export numeric_loop_matrix_and_vector!
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

export FiniteElements

export DiscreteModelWithFEMap
export GridWithFEMap
export add_mesh_displacement!
export update_coordinates!

export ConstantFESpace

include("FESpaceInterface.jl")

include("SingleFieldFESpaces.jl")

include("UnconstrainedFESpaces.jl")

include("ConformingFESpaces.jl")

include("Pullbacks.jl")

include("FESpaceFactories.jl")

include("PhysicalFEs.jl")

include("TrialFESpaces.jl")

include("Assemblers.jl")

include("SparseMatrixAssemblers.jl")

include("FEOperators.jl")

include("AffineFEOperators.jl")

include("FEOperatorsFromWeakForm.jl")

include("FESolvers.jl")

include("FEAutodiff.jl")

include("DiscontinuousFESpaces.jl")

include("FESpacesWithConstantFixed.jl")

include("ZeroMeanFESpaces.jl")

include("CLagrangianFESpaces.jl")

include("DirichletFESpaces.jl")

include("FESpacesWithLinearConstraints.jl")

include("DiscreteModelWithFEMaps.jl")

include("ConstantFESpaces.jl")

export get_free_values
"""
    get_free_values(args...)

!!! danger
    get_free_values has been removed. Use get_free_dof_values instead.
"""
function get_free_values(args...)
  @unreachable "get_free_values has been removed. Use get_free_dof_values instead."
end

export get_dirichlet_values
"""
    get_dirichlet_values(args...)

!!! danger
    get_dirichlet_values has been removed. Use get_dirichlet_dof_values instead.
"""
function get_dirichlet_values(args...)
  @unreachable "get_dirichlet_values has been removed. Use get_dirichlet_dof_values instead."
end

end # module
