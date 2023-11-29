"""
    struct ForwardEuler <: ODESolver end

Forward Euler ODE solver.
"""
struct ForwardEuler <: ODESolver
  sol::NonlinearSolver
  dt::Real
end

# ODESolver interface
function get_dt(solver::ForwardEuler)
  solver.dt
end

function allocate_dop_cache(
  solver::ForwardEuler,
  ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  nothing
end

function allocate_dop_cache(
  solver::ForwardEuler,
  ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
  r = allocate_residual(ode_op, t, (u, u), ode_cache)
  (J, r)
end

function allocate_sol_cache(solver::ForwardEuler)
  nothing
end

function DiscreteODEOperator(
  solver::ForwardEuler, ode_op::ODEOperator, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real
)
  ForwardEulerNonlinearOperator(ode_op, ode_cache, t0, u0, dt)
end

function DiscreteODEOperator(
  solver::ForwardEuler, ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real
)
  J, r = dop_cache
  ForwardEulerLinearOperator(ode_op, ode_cache, J, r, t0, u0, dt)
end

###############
# solve_step! #
###############
function solve_step!(
  uF::AbstractVector,
  solver::ForwardEuler, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real,
  cache
)
  # Unpack solver
  dt = get_dt(solver)

  # Allocate or unpack cache
  if isnothing(cache)
    ode_cache = allocate_cache(ode_op, t0, (u0, u0))
    dop_cache = allocate_dop_cache(solver, ode_op, ode_cache, t0, u0)
    sol_cache = allocate_sol_cache(solver)
  else
    ode_cache, dop_cache, sol_cache = cache
  end

  # Create discrete ODE operator
  dop = DiscreteODEOperator(
    solver, ode_op, ode_cache, dop_cache,
    t0, u0, dt
  )

  # Solve discrete ODE operator
  uF, sol_cache = solve_dop!(uF, dop, solver.sol, sol_cache)
  tF = t0 + dt

  # Update cache
  cache = (ode_cache, dop_cache, sol_cache)

  (uF, tF, cache)
end

"""
    struct ForwardEulerNonlinearOperator <: DiscreteODEOperator

Nonlinear discrete operator corresponding to the forward Euler scheme:
```math
residual(t_n, u_n, v) = 0.
```
"""
struct ForwardEulerNonlinearOperator <: DiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  t0::Real
  u0::AbstractVector
  dt::Real
end

function Algebra.allocate_residual(
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t0, op.u0
  allocate_residual(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t0, op.u0
  residual!(r, op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t0, op.u0
  allocate_jacobian(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t0, op.u0
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, op.ode_op, t, (u, v), 1, 1, op.ode_cache)
end

function solve_dop!(
  uF::AbstractVector,
  op::ForwardEulerNonlinearOperator, sol::NonlinearSolver, cache
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  t0, u0, dt = op.t0, op.u0, op.dt

  update_cache!(ode_cache, ode_op, t0)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  cache = solve!(vF, sol, op, cache)
  _u_from_v!(uF, u0, dt, vF)

  (uF, cache)
end

"""
    struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator

Linear discrete operator corresponding to the forward Euler scheme:
```math
residual(t_n, u_n, v) = mass(t_n, u_n) v + res(t_n, u_n) = 0.
```
"""
struct ForwardEulerLinearOperator <: LinearDiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  u0::AbstractVector
  dt::Real
end

Algebra.get_matrix(op::ForwardEulerLinearOperator) = op.J

Algebra.get_vector(op::ForwardEulerLinearOperator) = op.r

function solve_dop!(
  uF::AbstractVector,
  op::ForwardEulerLinearOperator, sol::NonlinearSolver, cache
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  J, r = op.J, op.r
  t0, u0, dt = op.t0, op.u0, op.dt

  update_cache!(ode_cache, ode_op, t0)

  u, v = u0, u0
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, ode_op, t0, (u, v), 1, 1, ode_cache)
  residual!(r, ode_op, t0, (u, v), ode_cache, include_highest=false)
  rmul!(r, -1)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  cache = solve!(vF, sol, op, cache)
  _u_from_v!(uF, u0, dt, vF)

  (uF, cache)
end
