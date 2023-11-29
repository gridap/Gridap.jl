"""
    struct ThetaMethod <: ODESolver end

θ-method ODE solver.
"""
struct ThetaMethod <: ODESolver
  sol::NonlinearSolver
  dt::Real
  θ::Real

  function ThetaMethod(sol, dt, θ)
    θ01 = clamp(θ, 0, 1)
    if θ01 != θ
      msg = """
      The parameter θ of the θ-method must lie between zero and one.
      Setting θ to $(θ01).
      """
      @warn msg
    end

    if iszero(θ01)
      ForwardEuler(sol, dt)
    else
      new(sol, dt, θ01)
    end
  end
end

MidPoint(sol, dt) = ThetaMethod(sol, dt, 0.5)
BackwardEuler(sol, dt) = ThetaMethod(sol, dt, 1)

# ODESolver interface
function get_dt(solver::ThetaMethod)
  solver.dt
end

function allocate_dop_cache(
  solver::ThetaMethod,
  ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  similar(u)
end

function allocate_dop_cache(
  solver::ThetaMethod,
  ode_op::ODEOperator{LinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
  r = allocate_residual(ode_op, t, (u, u), ode_cache)
  (J, r)
end

function allocate_sol_cache(solver::ThetaMethod)
  nothing
end

function DiscreteODEOperator(
  solver::ThetaMethod, ode_op::ODEOperator, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real, tθ::Real, dtθ::Real
)
  ulc = dop_cache
  ThetaMethodNonlinearOperator(ode_op, ode_cache, ulc, t0, u0, dt, tθ, dtθ)
end

function DiscreteODEOperator(
  solver::ThetaMethod, ode_op::ODEOperator{LinearODE}, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real, tθ::Real, dtθ::Real
)
  J, r = dop_cache
  ThetaMethodLinearOperator(ode_op, ode_cache, J, r, t0, u0, dt, tθ, dtθ)
end

###############
# solve_step! #
###############
function solve_step!(
  uF::AbstractVector,
  solver::ThetaMethod, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real,
  cache
)
  # Unpack solver
  dt = get_dt(solver)
  θ = solver.θ
  dtθ = θ * dt
  tθ = t0 + dtθ

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
    t0, u0, dt, tθ, dtθ
  )

  # Solve discrete ODE operator
  uF, sol_cache = solve_dop!(uF, dop, solver.sol, sol_cache)
  tF = t0 + dt

  # Update cache
  cache = (ode_cache, dop_cache, sol_cache)

  (uF, tF, cache)
end

"""
    struct ThetaMethodNonlinearOperator <: DiscreteODEOperator

Nonlinear discrete operator corresponding to the θ-method scheme:
```math
residual(t_θ, u_n + θ * dt * v, v) = 0.
```
"""
struct ThetaMethodNonlinearOperator <: DiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  ulc::AbstractVector
  t0::Real
  u0::AbstractVector
  dt::Real
  tθ::Real
  dtθ::Real
end

function Algebra.allocate_residual(
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  t, u = op.tθ, op.ulc
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  allocate_residual(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  t, u = op.tθ, op.ulc
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  residual!(r, op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  t, u = op.tθ, op.ulc
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  allocate_jacobian(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  t, u = op.tθ, op.ulc
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, op.ode_op, t, (u, v), (op.dtθ, 1), op.ode_cache)
end

function solve_dop!(
  uF::AbstractVector,
  op::ThetaMethodNonlinearOperator, sol::NonlinearSolver, cache
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  tθ, u0, dt = op.tθ, op.u0, op.dt

  update_cache!(ode_cache, ode_op, tθ)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  cache = solve!(vF, sol, op, cache)
  _u_from_v!(uF, u0, dt, vF)

  (uF, cache)
end

"""
    struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the θ-method scheme:
```math
residual(t_θ, u_n + θ * dt * v, v) = mass(t_θ) v + stiffness(t_θ) (u_n + θ * dt * v) + res(t_θ) = 0.
```
"""
struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  u0::AbstractVector
  dt::Real
  tθ::Real
  dtθ::Real
end

Algebra.get_matrix(op::ThetaMethodLinearOperator) = op.J

Algebra.get_vector(op::ThetaMethodLinearOperator) = op.r

function solve_dop!(
  uF::AbstractVector,
  op::ThetaMethodLinearOperator, sol::NonlinearSolver, cache
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  J, r = op.J, op.r
  tθ, u0, dtθ, dt = op.tθ, op.u0, op.dtθ, op.dt

  update_cache!(ode_cache, ode_op, tθ)

  u, v = u0, u0
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, ode_op, tθ, (u, v), (dtθ, 1), ode_cache)
  residual!(r, ode_op, tθ, (u, v), ode_cache, include_highest=false)
  rmul!(r, -1)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  cache = solve!(vF, sol, op, cache)
  _u_from_v!(uF, u0, dt, vF)

  (uF, cache)
end
