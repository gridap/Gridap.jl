"""
    struct ThetaMethod <: ODESolver end

θ-method ODE solver:
* θ = 0     Forward Euler
* θ = 1/2   Crank-Nicolson / MidPoint
* θ = 1     Backward Euler
"""
struct ThetaMethod <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  θ::Float64

  function ThetaMethod(nls, dt, θ)
    if θ <= 0
      ForwardEuler(nls, dt)
    elseif θ >= 1
      BackwardEulerEuler(nls, dt)
    else
      new(nls, dt, θ)
    end
  end
end

MidPoint(nls, dt) = ThetaMethod(nls, dt, 0.5)

function solve_step!(
  uF::AbstractVector,
  solver::ThetaMethod, op::ODEOperator,
  u0::AbstractVector, t0::Real,
  cache
)
  if isnothing(cache)
    ode_cache = allocate_cache(op)
    u̇θ = similar(u0)
    nls_cache = nothing
  else
    ode_cache, u̇θ, nls_cache = cache
  end

  dt, θ = solver.dt, solver.θ
  dtθ = θ * dt
  tθ = t0 + dtθ
  tF = t0 + dt

  # Update Dirichlet boundary conditions
  ode_cache = update_cache!(ode_cache, op, tθ)

  # Create and solve discrete ODE operator
  nl_op = ThetaMethodSolverOperator(op, ode_cache, tθ, dtθ, u0, u̇θ)
  nls_cache = solve!(uF, solver.nls, nl_op, nls_cache)

  # The operator above solves for uθ
  # Express uF in terms of uθ
  axpy!(θ - 1, u0, uF)
  rdiv!(uF, θ)

  # Update cache
  cache = (ode_cache, u̇θ, nls_cache)

  (uF, tF, cache)
end

"""
    struct ThetaMethodSolverOperator <: NonlinearOperator end

Nonlinear operator that represents the θ-method nonlinear operator at a
given time step, i.e., residual(t, u_n+θ, (u_n+θ - u_n) / dt)
"""
struct ThetaMethodSolverOperator <: NonlinearOperator
  ode_op::ODEOperator
  ode_cache
  tθ::Float64
  dtθ::Float64
  u0::AbstractVector
  u̇θ::AbstractVector
end

function Algebra.zero_initial_guess(op::ThetaMethodSolverOperator)
  u = similar(op.u0)
  fill!(u, zero(eltype(u)))
  u
end

function Algebra.allocate_residual(
  op::ThetaMethodSolverOperator,
  u::AbstractVector
)
  allocate_residual(op.ode_op, op.tθ, (u, u), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ThetaMethodSolverOperator,
  u::AbstractVector
)
  uθ = u
  u0, u̇θ, dtθ = op.u0, op.u̇θ, op.dtθ
  _discrete_time_derivative!(u̇θ, u0, uθ, dtθ)
  residual!(r, op.ode_op, op.tθ, (uθ, u̇θ), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ThetaMethodSolverOperator,
  u::AbstractVector
)
  allocate_jacobian(op.ode_op, op.tθ, (u, u), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ThetaMethodSolverOperator,
  u::AbstractVector
)
  fillstored!(J, zero(eltype(J)))

  uθ = u
  u0, u̇θ, dtθ = op.u0, op.u̇θ, op.dtθ
  _discrete_time_derivative!(u̇θ, u0, uθ, dtθ)
  jacobians!(J, op.ode_op, op.tθ, (uθ, u̇θ), (1, inv(dtθ)), op.ode_cache)
end
