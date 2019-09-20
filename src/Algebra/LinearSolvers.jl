module LinearSolvers

using Gridap.Helpers
using SparseArrays
using LinearAlgebra
using Test

export solve
export solve!
export LinearSolver
export SymbolicSetup
export NumericalSetup
export symbolic_setup
export numerical_setup
export numerical_setup!
export test_linear_solver
export LUSolver
export BackslashSolver

abstract type LinearSolver end

abstract type SymbolicSetup end

abstract type NumericalSetup end

symbolic_setup(::LinearSolver,mat::AbstractMatrix)::SymbolicSetup = @abstractmethod

numerical_setup(::SymbolicSetup,mat::AbstractMatrix)::NumericalSetup = @abstractmethod

numerical_setup!(::NumericalSetup,mat::AbstractMatrix) = @abstractmethod

solve!(x::AbstractVector,::NumericalSetup,b::AbstractVector) = @abstractmethod

function LinearSolver end

function solve(ls::LinearSolver,A::AbstractMatrix,b::AbstractVector)
  ss = symbolic_setup(ls,A)
  ns = numerical_setup(ss,A)
  x = similar(b)
  solve!(x,ns,b)
  x
end

function test_linear_solver(
  ls::LinearSolver,A::AbstractMatrix,b::AbstractVector,x::AbstractVector)

  y = solve(ls,A,b)
  @test x ≈ y

  ss = symbolic_setup(ls,A)
  ns = numerical_setup(ss,A)
  numerical_setup!(ns,A)
  solve!(y,ns,b)
  @test x ≈ y

end

"""
Wrapper of the LU solver available in julia
"""
struct LUSolver <: LinearSolver end

struct LUSymbolicSetup <: SymbolicSetup end

mutable struct LUNumericalSetup{F} <: NumericalSetup
  factors::F
end

symbolic_setup(::LUSolver,mat::AbstractMatrix) = LUSymbolicSetup()

numerical_setup(::LUSymbolicSetup,mat::AbstractMatrix) = LUNumericalSetup(lu(mat))

function numerical_setup!(ns::LUNumericalSetup, mat::AbstractMatrix)
  fac = lu(mat)
  ns.factors = fac
end

function solve!(
  x::AbstractVector,ns::LUNumericalSetup,b::AbstractVector)
  y = ns.factors\b # the allocation of y can be avoided
  x .= y
end

"""
Wrapper of the backslash solver available in julia
This is typically faster than LU for a single solve
"""
struct BackslashSolver <: LinearSolver end

struct BackslashSymbolicSetup <: SymbolicSetup end

mutable struct BackslashNumericalSetup{T<:AbstractMatrix} <: NumericalSetup
  A::T
end

symbolic_setup(::BackslashSolver,mat::AbstractMatrix) = BackslashSymbolicSetup()

numerical_setup(::BackslashSymbolicSetup,mat::AbstractMatrix) = BackslashNumericalSetup(mat)

function numerical_setup!(ns::BackslashNumericalSetup, mat::AbstractMatrix)
  ns.A = mat
end

function solve!(
  x::AbstractVector,ns::BackslashNumericalSetup,b::AbstractVector)
   copyto!(x, ns.A\b)
end

end
