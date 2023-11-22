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
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using Gridap.MultiField

const ε = 100 * eps()
const OneOrMoreOfType{T} = Union{T,Tuple{Vararg{T}}}
const OneOrMoreVectors = OneOrMoreOfType{AbstractVector}
const OneOrMoreFunctions = OneOrMoreOfType{Function}

include("TimeDerivatives.jl")

export time_derivative
export ∂t
export ∂tt

include("ODEOperators.jl")

export ODEOperatorType
export NonlinearODE
export MassLinearODE
export ConstantMassODE

export ODEOperator
export MassLinearODEOperator
export ConstantMassODEOperator

export get_order
export allocate_cache
export update_cache!
export allocate_residual
export residual!
export allocate_jacobian
export jacobian!
export jacobians!

export test_ode_operator

include("ODESolvers.jl")

export ODESolver
export solve_step!
export solve

export test_ode_solver

include("ODESolutions.jl")

export ODESolution
export GenericODESolution

export test_ode_solution

include("TransientFESpaces.jl")

export TransientTrialFESpace
export TransientMultiFieldFESpace

export test_transient_trial_fe_space

include("TransientCellFields.jl")

export TransientCellField

include("TransientFEOperators.jl")

export TransientFEOperator

export test_transient_fe_operator

include("TransientFESolutions.jl")

export TransientFESolution

export test_transient_fe_solution
export test_transient_fe_solver

# include("_DiffEqsWrappers.jl")

end # module ODEs

const GridapODEs = ODEs
