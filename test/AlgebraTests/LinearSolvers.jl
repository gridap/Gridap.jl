module LinearSolvers

using Gridap.Helpers
using SparseArrays
using LinearAlgebra
using Test

export LinearSolver
export SymbolicSetup
export NumericalSetup
export solve
export solve!
export test_linear_solver
export LUSolver

abstract type LinearSolver end

abstract type SymbolicSetup end

abstract type NumericalSetup end

SymbolicSetup(
  ::LinearSolver,mat::AbstractMatrix)::SymbolicSetup = @abstractmethod

NumericalSetup(
  ::SymbolicSetup,mat::AbstractMatrix)::NumericalSetup = @abstractmethod

solve!(
  x::AbstractVector,::NumericalSetup,A::AbstractMatrix,b::AbstractVector) = @abstractmethod

function solve(ls::LinearSolver,A::AbstractMatrix,b::AbstractVector)
  ss = SymbolicSetup(ls,A)
  ns = NumericalSetup(ss,A)
  x = similar(b)
  solve!(x,ns,A,b)
  x
end

function test_linear_solver(
  ls::LinearSolver,A::AbstractMatrix,b::AbstractVector,x::AbstractVector)

  y = solve(ls,A,b)
  @test x ≈ y

end


"""
Wrapper of the LU solver available in julia
"""
struct LUSolver <: LinearSolver end

struct LUSymbolicSetup <: SymbolicSetup end

struct LUNumericalSetup{F} <: NumericalSetup
  factors::F
end

SymbolicSetup(::LUSolver,mat::AbstractMatrix) = LUSymbolicSetup()

NumericalSetup(::LUSymbolicSetup,mat::AbstractMatrix) = LUNumericalSetup(lu(mat))

function solve!(
  x::AbstractVector,ns::LUNumericalSetup,A::AbstractMatrix,b::AbstractVector)
  y = ns.factors\b # the allocation of y can be avoided
  x .= y
end

end
