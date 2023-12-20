"""

The exported names are
$(EXPORTS)
"""
module TransientFETools

using Test
using DocStringExtensions

using Gridap.Helpers

export ∂t

import Gridap.ODEs.ODETools: ∂t, ∂tt
import Gridap.ODEs.ODETools: time_derivative

export TransientTrialFESpace
export TransientMultiFieldFESpace
export test_transient_trial_fe_space
import Gridap.Fields: evaluate
import Gridap.Fields: evaluate!
import Gridap.MultiField: MultiFieldFESpace
using Gridap.FESpaces: FESpace
using Gridap.FESpaces: SingleFieldFESpace
using Gridap.FESpaces: TrialFESpace
using Gridap.FESpaces: ZeroMeanFESpace
using Gridap.FESpaces: get_free_dof_values
using Gridap.FESpaces: get_dirichlet_dof_values
using Gridap.FESpaces: TrialFESpace!
using Gridap.FESpaces: HomogeneousTrialFESpace
using Gridap.FESpaces: jacobian

import Gridap.Geometry: Triangulation
import Gridap.CellData: Measure
using Gridap.FESpaces: ∫

export TransientFEOperator
export TransientAffineFEOperator
export TransientConstantFEOperator
export TransientConstantMatrixFEOperator
export TransientRungeKuttaFEOperator
export TransientIMEXRungeKuttaFEOperator
export TransientEXRungeKuttaFEOperator
using Gridap.FESpaces: Assembler
using Gridap.FESpaces: SparseMatrixAssembler
import Gridap.ODEs.ODETools: allocate_cache
import Gridap.ODEs.ODETools: update_cache!
import Gridap.ODEs.ODETools: ODEOperator
import Gridap.ODEs.ODETools: AffineODEOperator
import Gridap.ODEs.ODETools: ConstantODEOperator
import Gridap.ODEs.ODETools: ConstantMatrixODEOperator
import Gridap.ODEs.ODETools: allocate_residual
import Gridap.ODEs.ODETools: allocate_jacobian
import Gridap.ODEs.ODETools: residual!
import Gridap.ODEs.ODETools: jacobian!
import Gridap.ODEs.ODETools: jacobians!
import Gridap.ODEs.ODETools: lhs!
import Gridap.ODEs.ODETools: rhs!
import Gridap.ODEs.ODETools: explicit_rhs!
import Gridap.ODEs.ODETools: OperatorType
using Gridap.ODEs.ODETools: Nonlinear
using Gridap.ODEs.ODETools: Affine
using Gridap.ODEs.ODETools: Constant
using Gridap.ODEs.ODETools: ConstantMatrix
import Gridap.FESpaces: get_algebraic_operator
import Gridap.FESpaces: assemble_vector!
import Gridap.FESpaces: assemble_matrix_add!
import Gridap.FESpaces: allocate_vector
import Gridap.FESpaces: allocate_matrix
using Gridap.FESpaces: get_fe_basis
using Gridap.FESpaces: get_trial_fe_basis
using Gridap.FESpaces: collect_cell_vector
using Gridap.FESpaces: collect_cell_matrix
using Gridap.FESpaces: return_type
import Gridap.FESpaces: SparseMatrixAssembler
import Gridap.FESpaces: get_trial
import Gridap.FESpaces: get_test
using Gridap.ODEs.ODETools: test_ode_operator
export test_transient_fe_operator

import Gridap.FESpaces: FESolver
import Gridap.ODEs.ODETools: ODESolver
import Gridap.Algebra: solve
import Gridap.Algebra: solve!
import Gridap.ODEs.ODETools: solve_step!
export test_transient_fe_solver

export TransientFEFunction
import Gridap.FESpaces: FEFunction
import Gridap.FESpaces: SingleFieldFEFunction
import Gridap.FESpaces: EvaluationFunction
import Gridap.MultiField: MultiFieldFEFunction
import Gridap.MultiField: num_fields

export TransientFESolution
import Gridap.Algebra: solve
import Gridap.ODEs.ODETools: ODESolution
import Gridap.ODEs.ODETools: GenericODESolution
import Base: iterate
export test_transient_fe_solution

export TransientCellField
using Gridap.CellData: CellField
using Gridap.CellData: CellFieldAt
using Gridap.CellData: GenericCellField
using Gridap.MultiField: MultiFieldCellField
using Gridap.FESpaces: FEBasis
import Gridap.CellData: get_data
import Gridap.CellData: get_triangulation
import Gridap.CellData: DomainStyle
import Gridap.CellData: gradient
import Gridap.CellData: ∇∇
import Gridap.CellData: change_domain
import Gridap.FESpaces: BasisStyle
using Gridap.FESpaces: Constrained, UnConstrained, AssemblyStrategy
using Gridap.MultiField: ConsecutiveMultiFieldStyle, BlockSparseMatrixAssembler
import Gridap.MultiField: ConstraintStyle, MultiFieldStyle, BlockMultiFieldStyle
import Gridap.FESpaces: zero_free_values, has_constraints, SparseMatrixAssembler
import Gridap.FESpaces: get_dof_value_type, get_vector_type

using BlockArrays

include("TransientFESpaces.jl")

include("TransientCellField.jl")

include("TransientMultiFieldCellField.jl")

include("TransientFEOperators.jl")

include("ODEOperatorInterfaces.jl")

include("TransientFESolutions.jl")

# export FETerm
# function FETerm(args...)
#   Helpers.@unreachable """\n
#   Function FETerm has been removed. The API for specifying the weak form has changed significantly.
#   See the gridap/Tutorials repo for some examples of how to use the new API.
#   This error message will be deleted in future versions.
#   """
# end

end #module
