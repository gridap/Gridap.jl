#############
# ODESolver #
#############
"""
An ODE solver is a map that given (t_n, us_n) returns (t_n+1, us_n+1) and the
corresponding updated cache. Here `us = (u_n, ∂t(u_n), ..., ∂t^(N-1)(u_n))` is
a tuple containing the time derivatives up to the order of the ODE operator
minus one.
"""
abstract type ODESolver <: GridapType end

"""
Evolve the ODE operator for one time step using the ODE solver
"""
function solve_step!(
  usF::VecOrNTupleVec,
  solver::ODESolver, op::ODEOperator,
  us0::VecOrNTupleVec, t0::Real, cache
)
  @abstractmethod
end

function solve_step!(
  usF::VecOrNTupleVec,
  solver::ODESolver, op::ODEOperator,
  us0::VecOrNTupleVec, t0::Real
)
  solve_step!(usF, solver, op, us0, t0, nothing)
end

function solve(
  solver::ODESolver, op::ODEOperator,
  u0::T, t0::Real, tF::Real
) where {T}
  GenericODESolution{T}(solver, op, u0, t0, tF)
end

########
# Test #
########
function test_ode_solver(solver::ODESolver, op::ODEOperator, u0, t0, tf)
  solution = solve(solver, op, u0, t0, tf)
  test_ode_solution(solution)
end
