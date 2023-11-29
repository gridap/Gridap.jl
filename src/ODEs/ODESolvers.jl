#######################
# DiscreteODEOperator #
#######################
"""
    abstract type DiscreteODEOperator <: NonlinearOperator end

Discrete ODE operator corresponding to an `ODEOperator` and an `ODESolver`.

# Mandatory
- [`solve_dop!(uF, op, sol, cache)`](@ref)
"""
abstract type DiscreteODEOperator <: NonlinearOperator end

"""
    solve_dop!(
      usF::OneOrMoreVectors,
      op::DiscreteODEOperator, sol::NonlinearSolver, cache
    ) -> (OneOrMoreVectors, CacheType)

Solve the discrete ODE operator.
"""
function solve_dop!(
  usF::OneOrMoreVectors,
  op::DiscreteODEOperator, sol::NonlinearSolver, cache
)
  @abstractmethod
end

#############################
# LinearDiscreteODEOperator #
#############################
"""
    abstract type LinearDiscreteODEOperator <: DiscreteODEOperator end

Discrete linear ODE operator corresponding to an `ODEOperator` and an
`ODESolver`.

# Mandatory
- [`get_matrix(op)`](@ref)
- [`get_vector(op)`](@ref)
"""
abstract type LinearDiscreteODEOperator <: DiscreteODEOperator end

function Algebra.get_matrix(op::LinearDiscreteODEOperator)
  @abstractmethod
end

function Algebra.get_vector(op::LinearDiscreteODEOperator)
  @abstractmethod
end

function Algebra.allocate_residual(
  op::LinearDiscreteODEOperator, v::AbstractVector
)
  similar(v)
end

function Algebra.residual!(
  r::AbstractVector,
  op::LinearDiscreteODEOperator, v::AbstractVector
)
  mul!(r, get_matrix(op), v)
  r .-= get_vector(op)
  r
end

function Algebra.allocate_jacobian(
  op::LinearDiscreteODEOperator, v::AbstractVector
)
  get_matrix(op)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::LinearDiscreteODEOperator, v::AbstractVector
)
  copy_entries!(J, get_matrix(op))
  J
end

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
- [`get_dt(solver)`](@ref)
- [`allocate_dop_cache(solver, ode_op, ode_cache, t, u)`](@ref)
- [`allocate_sol_cache(solver)`](@ref)
- [`DiscreteODEOperator(solver, ode_op, ode_cache, dop_cache, args...)`](@ref)
- [`solve_step!(usF, solver, op, us0, t0[, cache])`](@ref)

# Optional
- [`solve(solver, op, us0, t0, tF)`](@ref)
"""
abstract type ODESolver <: GridapType end

"""
    get_dt(solver::ODESolver) -> Real

Return the time step of the `ODESolver`.
"""
function get_dt(solver::ODESolver)
  @abstractmethod
end

"""
    allocate_dop_cache(
      solver::ODESolver,
      ode_op::ODEOperator, ode_cache,
      t::Real, u::AbstractVector
    ) -> CacheType

Allocate the cache of the discrete disctem of the `ODESolver`.
"""
function allocate_dop_cache(
  solver::ODESolver,
  ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  @abstractmethod
end

"""
    allocate_sol_cache(solver::ODESolver) -> CacheType

Allocate the cache of the solver for the discrete disctem of the `ODESolver`.
"""
function allocate_sol_cache(solver::ODESolver)
  @abstractmethod
end

"""
    DiscreteODEOperator(
      solver::ODESolver, ode_op::ODEOperator, ode_cache,
      dop_cache, args...
    ) -> DiscreteODEOperator

Return the discrete ODE operator corresponding to the `ODEOperator` and the
`ODESolver`.
"""
function DiscreteODEOperator(
  solver::ODESolver, ode_op::ODEOperator, ode_cache,
  dop_cache, args...
)
  @abstractmethod
end

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
  solver::ODESolver, ode_op::ODEOperator,
  t0::Real, us0::OneOrMoreVectors, args...
)
  @test get_dt(solver) isa Real

  u0 = (us0 isa Tuple) ? first(us0) : us0
  ode_cache = allocate_cache(ode_op, t0, us0)
  dop_cache = allocate_dop_cache(solver, ode_op, ode_cache, t0, u0)

  dop = DiscreteODEOperator(solver, ode_op, ode_cache, dop_cache, args...)
  @test dop isa DiscreteODEOperator

  usF = copy.(us0)
  (usF, tF, cache) = solve_step!(usF, solver, ode_op, us0, t0)
  (usF, tF, cache) = solve_step!(usF, solver, ode_op, us0, t0, cache)

  @test usF isa OneOrMoreVectors
  @test tF isa Real

  true
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

include("ODESolvers/RungeKutta.jl")

include("ODESolvers/GeneralizedAlpha.jl")
