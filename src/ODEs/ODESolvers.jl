#############
# ODESolver #
#############
"""
    abstract type ODESolver <: GridapType end

An ODE solver is a map that given (t_n, us_n) returns (t_n+1, us_n+1) and the
corresponding updated cache. Here `us = (u_n, ∂t(u_n), ..., ∂t^(N-1)(u_n))` is
a tuple containing the time derivatives up to the order of the ODE operator
minus one.

- [`solve_step!(usF, solver, op, us0, t0[, cache])`](@ref)
- [`solve(solver, op, us0, t0, tF)`](@ref)
"""
abstract type ODESolver <: GridapType end

"""
    solve_step!(
      usF::OneOrMoreVectors,
      solver::ODESolver, op::ODEOperator,
      us0::OneOrMoreVectors, t0::Real
      [, cache]
    ) -> Tuple{Real,OneOrMoreVectors,CacheType}

Perform one time step of the ODE operator from `t0` with initial state `us0`
with the ODE solver.
"""
function solve_step!(
  usF::OneOrMoreVectors,
  solver::ODESolver, op::ODEOperator,
  us0::OneOrMoreVectors, t0::Real,
  cache
)
  @abstractmethod
end

function solve_step!(
  usF::OneOrMoreVectors,
  solver::ODESolver, op::ODEOperator,
  us0::OneOrMoreVectors, t0::Real
)
  solve_step!(usF, solver, op, us0, t0, nothing)
end

"""
    solve(
      solver::ODESolver, op::ODEOperator,
      us0::OneOrMoreVectors, t0::Real, tF::Real
    ) -> ODESolution

Create an `ODESolution` wrapper around the ODE operator and ODE solver,
starting with state `us0` at time `t0`, to be evolved until `tF`
"""
function Algebra.solve(
  solver::ODESolver, op::ODEOperator,
  us0::OneOrMoreVectors, t0::Real, tF::Real
)
  T = typeof(us0)
  GenericODESolution{T}(solver, op, us0, t0, tF)
end

########
# Test #
########
"""
    test_ode_solver(
      solver::ODESolver, op::ODEOperator,
      us0::OneOrMoreVectors, t0::Real, tf::Real
    ) -> Bool

Test the interface of `ODESolver` specializations
"""
function test_ode_solver(
  solver::ODESolver, op::ODEOperator,
  us0::OneOrMoreVectors, t0::Real, tf::Real
)
  solution = solve(solver, op, us0, t0, tf)
  test_ode_solution(solution)
end

##################
# Import solvers #
##################
function _discrete_time_derivative!(u̇F, u0, u, dt)
  copy!(u̇F, u)
  axpy!(-1, u0, u̇F)
  rdiv!(u̇F, dt)
  u̇F
end

include("ODESolvers/ForwardEuler.jl")

include("ODESolvers/BackwardEuler.jl")

include("ODESolvers/ThetaMethod.jl")

include("ODESolvers/Tableaus.jl")
