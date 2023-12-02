"""

The exported names are
$(EXPORTS)
"""
module ODEs

using Test
using DocStringExtensions

using LinearAlgebra
using LinearAlgebra: fillstored!
using ForwardDiff
using BlockArrays

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.CellData: OperationCellField
using Gridap.FESpaces
using Gridap.FESpaces: SingleFieldFEBasis
using Gridap.MultiField

const ε = 100 * eps()

include("TimeDerivatives.jl")

export time_derivative
export ∂t
export ∂tt

include("ODEOperators.jl")

export ODEOperatorType
export NonlinearODE
export MassLinearODE
export LinearODE

export ODEOperator
export MassLinearODEOperator
export LinearODEOperator

export allocate_odeopcache
export update_odeopcache!
export jacobians!
export is_jacobian_constant
export is_forcing_constant

export test_ode_operator

include("ODESolvers.jl")

export DiscreteODEOperator
export LinearDiscreteODEOperator

export ODESolver
export get_dt
export allocate_disopcache
export allocate_disslvrcache
export solve_step!

export test_ode_solver

export ForwardEuler

export ThetaMethod
export MidPoint
export BackwardEuler

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

export TableauName
export ButcherTableau
export available_tableaus

export RungeKutta

export GeneralizedAlpha

include("ODESolutions.jl")

export ODESolution
export GenericODESolution

export test_ode_solution

include("TransientFESpaces.jl")

export AbstractTransientTrialFESpace
export allocate_space

export TransientTrialFESpace
export TransientMultiFieldFESpace

export test_transient_trial_fe_space

include("TransientCellFields.jl")

export TransientCellField
export TransientSingleFieldCellField
export TransientMultiFieldCellField
export TransientFEBasis

include("TransientFEOperators.jl")

export TransientFEOperator
export allocate_feopcache
export update_feopcache!
export get_assembler
export get_res
export get_jacs
export get_mass
export get_forms

export TransientFEOpFromWeakForm
export TransientMassLinearFEOpFromWeakForm
export TransientMassLinearFEOperator
export TransientLinearFEOpFromWeakForm
export TransientLinearFEOperator

export test_transient_fe_operator

include("ODEOperatorsFromFEOperators.jl")

export ODEOpFromFEOpCache
export ODEOpFromFEOp

include("TransientFESolutions.jl")

export TransientFESolution

export test_transient_fe_solution
export test_transient_fe_solver

# include("_DiffEqsWrappers.jl")

end # module ODEs

# TODO useful?
const GridapODEs = ODEs
