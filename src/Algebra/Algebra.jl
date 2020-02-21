"""

The exported names are
$(EXPORTS)
"""
module Algebra

using NLsolve
using DocStringExtensions
using SparseArrays
using LinearAlgebra
using Test

using Gridap.Helpers

export push_coo!
export is_entry_stored
export finalize_coo!
export sparse_from_coo
export add_entry!
export fill_entries!
export copy_entries!
export create_coo_vectors
export allocate_coo_vectors
export sparsecsr
export symsparsecsr
export hasrowmajororder
export hascolmajororder
export colvals
export getptr
export getindices
export SparseMatrixCSR
export SymSparseMatrixCSR

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
export allocate_residual
export allocate_jacobian
export allocate_residual_and_jacobian
export zero_initial_guess
export test_non_linear_operator

export NonLinearSolver
export test_non_linear_solver

export NewtonRaphsonSolver
export NLSolver

export AffineOperator
export get_matrix
export get_vector

export SparseMatrixCSR
export SymSparseMatrixCSR

include("SparseMatrices.jl")

include("CompressedSparseMatrices.jl")

include("NonLinearOperators.jl")

include("NonLinearSolvers.jl")

include("LinearSolvers.jl")

include("NLSolvers.jl")

end # module
