"""
    struct ForwardEuler <: ODESolver end

Forward Euler ODE solver
"""
struct ForwardEuler <: ODESolver
  nls::NonlinearSolver
  dt::Float64
end

function solve_step!(
  uF::AbstractVector,
  solver::ForwardEuler, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real,
  cache
)

  if isnothing(cache)
    ode_cache = allocate_cache(ode_op)
    u̇F = similar(u0)
    nls_cache = nothing
  else
    ode_cache, u̇F, nls_cache = cache
  end

  dt = solver.dt
  tF = t0 + dt

  # Update Dirichlet boundary conditions
  ode_cache = update_cache!(ode_cache, ode_op, t0)

  # Create and solve discrete ODE operator
  nl_op = ForwardEulerSolverOperator(ode_op, ode_cache, t0, dt, u0, u̇F)
  nls_cache = solve!(uF, solver.nls, nl_op, nls_cache)

  # Update cache
  cache = (ode_cache, u̇F, nls_cache)

  (uF, tF, cache)
end

"""
    struct ForwardEulerSolverOperator <: NonlinearOperator

Nonlinear operator that represents the Forward Euler nonlinear operator at a
given time step, i.e., residual(t, u_n, (u_n+1 - u_n) / dt)
"""
struct ForwardEulerSolverOperator <: NonlinearOperator
  ode_op::ODEOperator
  ode_cache
  t0::Float64
  dt::Float64
  u0::AbstractVector
  u̇F::AbstractVector
end

function Algebra.zero_initial_guess(op::ForwardEulerSolverOperator)
  u = similar(op.u0)
  fill!(u, zero(eltype(u)))
  u
end

function Algebra.allocate_residual(
  op::ForwardEulerSolverOperator,
  u::AbstractVector
)
  allocate_residual(op.ode_op, op.t0, (u, u), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ForwardEulerSolverOperator,
  u::AbstractVector
)
  u0, u̇F, dt = op.u0, op.u̇F, op.dt
  _discrete_time_derivative!(u̇F, u0, u, dt)
  residual!(r, op.ode_op, op.t0, (u0, u̇F), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ForwardEulerSolverOperator,
  u::AbstractVector
)
  allocate_jacobian(op.ode_op, op.t0, (u, u), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ForwardEulerSolverOperator,
  u::AbstractVector
)
  fillstored!(J, zero(eltype(J)))

  u0, u̇F, dt = op.u0, op.u̇F, op.dt
  _discrete_time_derivative!(u̇F, u0, u, dt)
  jacobians!(J, op.ode_op, op.t0, (u0, u̇F), (0, 1 / dt), op.ode_cache)
end
