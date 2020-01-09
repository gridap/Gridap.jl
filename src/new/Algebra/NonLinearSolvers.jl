
"""
    abstract type NonLinearSolver <: GridapType end

- [`solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator)`](@ref)
- [`solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator, cache)`](@ref)
"""
abstract type NonLinearSolver <: GridapType end

"""
    solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator)

Usage:

   cache = solve!(x,nls,op)

The returned `cache` object can be used in subsequent solves:

   solve!(x,nls,op,cache)
"""
function solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator)
  @abstractmethod
end

"""
    solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator,cache)

Solve using the cache object from a previous solve.
"""
function solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator,cache)
  @abstractmethod
end


"""
    solve(nls::NonLinearSolver,op::NonLinearOperator)

Creates and uses a zero initial guess.
"""
function solve(nls::NonLinearSolver,op::NonLinearOperator)
  x = allocate_solution(op)
  fill!(x,zero(eltype(x)))
  solve!(x,nls,op)
  x
end

"""
    test_non_linear_solver(
      nls::NonLinearSolver,
      op::NonLinearOperator,
      x0::AbstractVector,
      x::AbstractVector,
      pred::Function=isapprox)
"""
function test_non_linear_solver(
  nls::NonLinearSolver,
  op::NonLinearOperator,
  x0::AbstractVector,
  x::AbstractVector,
  pred::Function=isapprox)

  x1 = copy(x0)
  cache = solve!(x1,nls,op)
  @test pred(x1,x)

  x1 = copy(x0)
  solve!(x1,nls,op,cache)
  @test pred(x1,x)

end


