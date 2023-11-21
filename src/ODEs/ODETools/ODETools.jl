"""

The exported names are
$(EXPORTS)
"""
module ODETools

using Test
using DocStringExtensions

import Base: iterate
using ForwardDiff

using Gridap.Helpers: @abstractmethod, GridapType

import Gridap.Algebra: allocate_residual
import Gridap.Algebra: residual!
import Gridap.Algebra: allocate_jacobian
import Gridap.Algebra: jacobian!

using Gridap.Arrays: get_array

using Gridap.TensorValues: VectorValue, TensorValue

using Gridap.Fields: return_type

const ϵ = 100 * eps()
const VecOrNTupleVec = Union{AbstractVector,Tuple{Vararg{AbstractVector}}}

include("DiffOperators.jl")

export ∂t
export ∂tt
export time_derivative

include("ODEOperators.jl")

export ODEOperatorType
export NonlinearODE
export MassLinearODE
export ConstantMassODE
export ODEOperator
export MassLinearODEOperator
export ConstantMassODEOperator
export get_order
export jacobians!
export allocate_cache
export update_cache!
export test_ode_operator

include("ODESolvers.jl")

export ODESolver
export solve_step!
export test_ode_solver

include("ODESolutions.jl")

export ODESolution
export GenericODESolution
export test_ode_solution

end # module ODETools
