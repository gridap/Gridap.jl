"""
    struct ThetaMethod <: ODESolver end

θ-method ODE solver.
"""
struct ThetaMethod <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  θ::Float64

  function ThetaMethod(nls, dt, θ)
    θ01 = clamp(θ, 0, 1)
    if θ01 != θ
      msg = """
      The parameter θ of ThetaMethod must be between zero and one.
      Setting θ to $(θ01).
      """
      println(msg)
    end

    if iszero(θ01)
      ForwardEuler(nls, dt)
    else
      new(nls, dt, θ01)
    end
  end
end

MidPoint(nls, dt) = ThetaMethod(nls, dt, 0.5)
BackwardEuler(nls, dt) = ThetaMethod(nls, dt, 1)

function solve_step!(
  uF::AbstractVector,
  solver::ThetaMethod, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real,
  cache
)
  if isnothing(cache)
    ode_cache = allocate_cache(ode_op, t0, (u0, u0))
    nls_cache = nothing
    solver_cache = allocate_solver_cache(solver, ode_op, ode_cache, t0, u0)
  else
    ode_cache, nls_cache, solver_cache = cache
  end

  dt, θ = solver.dt, solver.θ
  dtθ = θ * dt
  tθ = t0 + dtθ
  tF = t0 + dt

  # Update Dirichlet boundary conditions
  ode_cache = update_cache!(ode_cache, ode_op, tθ)

  # Create and solve discrete ODE operator
  nl_op = ThetaMethodOperator(
    ode_op, ode_cache, solver_cache,
    tθ, u0, dtθ
  )
  v = uF
  nls_cache = solve!(v, solver.nls, nl_op, nls_cache)
  uF = _u_from_v!(uF, u0, dt, v)

  # Update cache
  cache = (ode_cache, nls_cache, solver_cache)

  (uF, tF, cache)
end

function allocate_solver_cache(
  solver::ThetaMethod, ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  (similar(u),)
end

function allocate_solver_cache(
  solver::ThetaMethod, ode_op::ODEOperator{LinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
  r = allocate_residual(ode_op, t, (u, u), ode_cache)
  (J, r)
end

"""
    ThetaMethodOperator(
      ode_op::ODEOperator, ode_cache, solver_cache,
      tθ::Real, u0::AbstractVector, dtθ::Real
    ) -> Union{ThetaMethodNonlinearOperator, ThetaMethodLinearOperator}

Return the linear or nonlinear Theta Method operator, depending on the
type of `ODEOperator`.
"""
function ThetaMethodOperator(
  ode_op::ODEOperator, ode_cache, solver_cache,
  tθ::Real, u0::AbstractVector, dtθ::Real
)
  u, = solver_cache
  ThetaMethodNonlinearOperator(
    ode_op, ode_cache,
    tθ, u0, dtθ, u
  )
end

function ThetaMethodOperator(
  ode_op::ODEOperator{LinearODE}, ode_cache, solver_cache,
  tθ::Real, u0::AbstractVector, dtθ::Real
)
  ThetaMethodLinearOperator(
    ode_op, ode_cache, solver_cache,
    tθ, u0, dtθ
  )
end

"""
    struct ThetaMethodNonlinearOperator <: NonlinearOperator

Theta method nonlinear operator at a given time step, i.e.
```math
residual(t_n+θ, u_n+θ, (u_n+θ - u_n) / θ / dt) = 0.
```
"""
struct ThetaMethodNonlinearOperator <: NonlinearOperator
  ode_op::ODEOperator
  ode_cache
  tθ::Float64
  u0::AbstractVector
  dtθ::Float64
  u::AbstractVector
end

function Algebra.zero_initial_guess(op::ThetaMethodNonlinearOperator)
  v = similar(op.u0)
  fill!(v, zero(eltype(v)))
  v
end

function Algebra.allocate_residual(
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  u = op.u
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  allocate_residual(op.ode_op, op.tθ, (u, v), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  u = op.u
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  residual!(r, op.ode_op, op.tθ, (u, v), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  u = op.u
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  allocate_jacobian(op.ode_op, op.tθ, (u, v), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ThetaMethodNonlinearOperator,
  v::AbstractVector
)
  fillstored!(J, zero(eltype(J)))
  u = op.u
  u = _u_from_v!(u, op.u0, op.dtθ, v)
  jacobians!(J, op.ode_op, op.tθ, (u, v), (op.dtθ, 1), op.ode_cache)
end

"""
    ThetaMethodLinearOperator(
      ode_op, ode_cache, solver_cache,
      tθ, u0, dtθ
    ) -> AffineOperator

Theta method linear operator at a given time step, i.e.
```math
mass(t_n+θ) * (u_n+θ - u_n) / θ / dt + stiffness(t_n+θ) * u_n+θ + res(t_n+θ) = 0.
```
"""
function ThetaMethodLinearOperator(
  ode_op, ode_cache, solver_cache,
  tθ, u0, dtθ
)
  J, r = solver_cache

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, ode_op, tθ, (u0, u0), (dtθ, 1), ode_cache)
  residual!(r, ode_op, tθ, (u0, u0), ode_cache, include_highest=false)
  rmul!(r, -1)

  AffineOperator(J, r)
end
