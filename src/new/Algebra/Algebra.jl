"""

The exported names are
$(EXPORTS)
"""
module Algebra

using DocStringExtensions
using SparseArrays
using LinearAlgebra
using Test

using Gridap.Helpers

export LinearSolver
export SymbolicSetup
export NumericalSetup
export solve
export solve!
export symbolic_setup
export numerical_setup
export numerical_setup!
export test_linear_solver

export LUSolver
export BackslashSolver

export NonLinearOperator
export residual!
export residual
export jacobian!
export jacobian
export residual_and_jacobian!
export residual_and_jacobian
export num_domain_dims
export num_range_dims
export allocate_solution
export allocate_residual
export allocate_jacobian
export allocate_residual_and_jacobian
export zero_initial_guess
export test_non_linear_operator
export is_square

export NonLinearSolver
export test_non_linear_solver

export NewtonRaphsonSolver

include("LinearSolvers.jl")

include("NonLinearOperators.jl")

include("NonLinearSolvers.jl")

end # module
