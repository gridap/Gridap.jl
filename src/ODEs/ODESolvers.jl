#############
# ODESolver #
#############
"""
    abstract type ODESolver <: GridapType end

An `ODESolver` is a map that given (t_n, us_n) returns (t_n+1, us_n+1) and the
corresponding updated cache. Here `us_n` is a vector of size `N-1`, where `N` is
the order of the `ODEOperator`, and `us_n[k] = ∂t^k(u)(t_n)` is the `k`-th-order
time derivative of `u` at `t_n`.

# Mandatory
- [`solve_step!(usF, solver, op, us0, t0[, cache])`](@ref)

# Optional
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

Perform one time step of the `ODEOperator` with the `ODESolver` from `t0` with
initial state `us0`.
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

Create an `ODESolution` wrapper around the `ODEOperator` and `ODESolver`,
starting with state `us0` at time `t0`, to be evolved until `tF`.
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
      us0::OneOrMoreVectors, t0::Real, tF::Real
    ) -> Bool

Test the interface of `ODESolver` specializations.
"""
function test_ode_solver(
  solver::ODESolver, op::ODEOperator,
  us0::OneOrMoreVectors, t0::Real, tF::Real
)
  solution = solve(solver, op, us0, t0, tF)
  test_ode_solution(solution)
end

##################
# Import solvers #
##################
"""
    _v_from_u(
      v::AbstractVector,
      u::AbstractVector, u0::AbstractVector, dt::Real
    ) -> AbstractVector

Safely write `(u - u0) / dt` into `v`.
"""
function _v_from_u(
  v::AbstractVector,
  u::AbstractVector, u0::AbstractVector, dt::Real
)
  if u !== v
    copy!(v, u)
  end
  axpy!(-1, u0, v)
  rdiv!(v, dt)
  v
end

"""
    _u_from_v!(
      u::AbstractVector,
      u0::AbstractVector, dt::Real, v::AbstractVector
    ) -> AbstractVector

Safely write `u0 + dt * v` into `u`.
"""
function _u_from_v!(
  u::AbstractVector,
  u0::AbstractVector, dt::Real, v::AbstractVector
)
  if u === v
    axpby!(1, u0, dt, u)
  else
    copy!(u, u0)
    axpy!(dt, v, u)
  end
  u
end

include("ODESolvers/ForwardEuler.jl")

include("ODESolvers/ThetaMethod.jl")

include("ODESolvers/Tableaus.jl")
