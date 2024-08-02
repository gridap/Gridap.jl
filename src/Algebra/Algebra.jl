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
using SparseMatricesCSR
using BlockArrays

using Gridap.Helpers

import Base: convert, size, getindex, show, count, *
import LinearAlgebra: mul!
import SparseArrays: nnz, nonzeros, nzrange, findnz, rowvals

export length_to_ptrs!
export rewind_ptrs!
export allocate_vector
export allocate_matrix
export allocate_matrix_and_vector
export allocate_in_domain
export allocate_in_range
export add_entries!
export muladd!
export axpy_entries!
export nz_counter
export nz_allocation
export create_from_nz
export LoopStyle
export Loop
export DoNotLoop
export ArrayBuilder
export SparseMatrixBuilder
export MinMemory
export MinCPU
export get_array_type

export nz_index
export is_entry_stored
export finalize_coo!
export sparse_from_coo
export add_entry!
export copy_entries!
export allocate_coo_vectors
export push_coo!

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

export NonlinearOperator
export residual!
export residual
export jacobian!
export jacobian
export hessian
export hessian!
export residual_and_jacobian!
export residual_and_jacobian
export allocate_residual
export allocate_jacobian
export allocate_residual_and_jacobian
export zero_initial_guess
export test_nonlinear_operator

export NonlinearSolver
export test_nonlinear_solver

export NewtonRaphsonSolver
export NLSolver

export AffineOperator
export get_matrix
export get_vector

include("AlgebraInterfaces.jl")

include("SparseMatrixCSC.jl")

include("SparseMatrixCSR.jl")

include("SymSparseMatrixCSR.jl")

include("NonlinearOperators.jl")

include("NonlinearSolvers.jl")

include("LinearSolvers.jl")

include("NLSolvers.jl")

end # module
