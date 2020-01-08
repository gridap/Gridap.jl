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

include("LinearSolvers.jl")


end # module
