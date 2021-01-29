
"""
    abstract type NonlinearSolver <: GridapType end

- [`solve!(x::AbstractVector,nls::NonlinearSolver,op::NonlinearOperator)`](@ref)
- [`solve!(x::AbstractVector,nls::NonlinearSolver,op::NonlinearOperator, cache)`](@ref)
"""
abstract type NonlinearSolver <: GridapType end

"""
    solve!(x::AbstractVector,nls::NonlinearSolver,op::NonlinearOperator)

Usage:

   cache = solve!(x,nls,op)

The returned `cache` object can be used in subsequent solves:

   cache = solve!(x,nls,op,cache)
"""
function solve!(x::AbstractVector,nls::NonlinearSolver,op::NonlinearOperator)
  solve!(x,nls,op,nothing)
end
function solve!(x::Tuple{AbstractVector,AbstractVector},nls::NonlinearSolver,op::NonlinearOperator)
  solve!(x,nls,op,nothing)
end

"""
    solve!(x::AbstractVector,nls::NonlinearSolver,op::NonlinearOperator,cache)

Solve using the cache object from a previous solve.
"""
function solve!(x::AbstractVector,nls::NonlinearSolver,op::NonlinearOperator,cache)
  @abstractmethod
end


"""
    solve(nls::NonlinearSolver,op::NonlinearOperator)

Creates and uses a zero initial guess.
"""
function solve(nls::NonlinearSolver,op::NonlinearOperator)
  x = zero_initial_guess(op)
  solve!(x,nls,op)
  x
end

"""
    test_nonlinear_solver(
      nls::NonlinearSolver,
      op::NonlinearOperator,
      x0::AbstractVector,
      x::AbstractVector,
      pred::Function=isapprox)
"""
function test_nonlinear_solver(
  nls::NonlinearSolver,
  op::NonlinearOperator,
  x0::AbstractVector,
  x::AbstractVector,
  pred::Function=isapprox)

  x1 = copy(x0)
  cache = solve!(x1,nls,op)
  @test pred(x1,x)

  x1 = copy(x0)
  cache = solve!(x1,nls,op,cache)
  @test pred(x1,x)

  x1 = copy(x0)
  cache = solve!(x1,nls,op,cache)
  @test pred(x1,x)

end
