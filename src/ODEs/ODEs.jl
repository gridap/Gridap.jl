"""

The exported names are
$(EXPORTS)
"""
module ODEs

using Test
using DocStringExtensions

using LinearAlgebra
using LinearAlgebra: fillstored!
using SparseArrays
using BlockArrays
using NLsolve
using ForwardDiff

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Algebra: NLSolversCache
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.CellData: OperationCellField
using Gridap.CellData: CellFieldAt
using Gridap.FESpaces
using Gridap.FESpaces: SingleFieldFEBasis
using Gridap.MultiField

const ε = 100 * eps()

include("TimeDerivatives.jl")

export TimeSpaceFunction
export time_derivative
export ∂t
export ∂tt

include("ODEOperators.jl")

export ODEOperatorType
export NonlinearODE
export AbstractQuasilinearODE
export QuasilinearODE
export AbstractSemilinearODE
export SemilinearODE
export AbstractLinearODE
export LinearODE

export ODEOperator
export get_num_forms
export get_forms
export is_form_constant
export allocate_odeopcache
export update_odeopcache!
export jacobian_add!

export IMEXODEOperator
export get_imex_operators

export GenericIMEXODEOperator

export test_ode_operator

include("StageOperators.jl")

export StageOperator
export NonlinearStageOperator
export LinearStageOperator

export massless_residual_weights

include("ODESolvers.jl")

export ODESolver
export allocate_odecache
export ode_start
export ode_march!
export ode_finish!

export test_ode_solver

export ForwardEuler

export ThetaMethod
export MidPoint
export BackwardEuler

export GeneralizedAlpha1

export TableauType
export ExplicitTableau
export ImplicitTableau
export FullyImplicitTableau
export DiagonallyImplicitTableau
export ImplicitExplicitTableau

export AbstractTableau
export GenericTableau
export EmbeddedTableau
export get_embedded_weights
export get_embedded_order
export IMEXTableau
export get_imex_tableaus
export is_padded

export TableauName
export ButcherTableau
export available_tableaus
export available_imex_tableaus

export RungeKutta

export GeneralizedAlpha2
export Newmark

include("ODESolutions.jl")

export ODESolution
export GenericODESolution

export test_ode_solution

include("TransientFESpaces.jl")

export allocate_space

export TransientTrialFESpace
export TransientMultiFieldFESpace

export test_tfe_space

include("TransientCellFields.jl")

export TransientCellField
export TransientSingleFieldCellField
export TransientMultiFieldCellField
export TransientFEBasis

include("TransientFEOperators.jl")

export TransientFEOperator
export get_assembler
export get_res
export get_jacs
export allocate_tfeopcache
export update_tfeopcache!

export TransientFEOpFromWeakForm
export TransientQuasilinearFEOpFromWeakForm
export TransientQuasilinearFEOperator
export TransientSemilinearFEOpFromWeakForm
export TransientSemilinearFEOperator
export TransientLinearFEOpFromWeakForm
export TransientLinearFEOperator

export TransientIMEXFEOperator
export GenericTransientIMEXFEOperator

export test_tfe_operator

include("ODEOpsFromTFEOps.jl")

export ODEOpFromTFEOpCache
export ODEOpFromTFEOp

include("TransientFESolutions.jl")

export TransientFESolution

export test_tfe_solution
export test_tfe_solver

# include("_DiffEqsWrappers.jl")

end # module ODEs

# TODO useful?
const GridapODEs = ODEs
