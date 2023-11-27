"""
    struct BackwardEuler <: ODESolver end

BackwardEuler ODE solver
"""
struct BackwardEuler <: ODESolver
  nls::NonlinearSolver
  dt::Float64
end

function solve_step!(
  uF::AbstractVector,
  solver::BackwardEuler, op::ODEOperator,
  u0::AbstractVector, t0::Real,
  cache
)
  if isnothing(cache)
    ode_cache = allocate_cache(op)
    u̇F = similar(u0)
    nls_cache = nothing
  else
    ode_cache, u̇F, nls_cache = cache
  end

  dt = solver.dt
  tF = t0 + dt

  # Update Dirichlet boundary conditions
  ode_cache = update_cache!(ode_cache, op, tF)

  # Create and solve discrete ODE operator
  nl_op = BackwardEulerSolverOperator(op, ode_cache, tF, dt, u0, u̇F)
  nls_cache = solve!(uF, solver.nls, nl_op, nls_cache)

  # Update cache
  cache = (ode_cache, u̇F, nls_cache)

  (uF, tF, cache)
end

"""
    struct BackwardEulerSolverOperator <: NonlinearOperator end

Nonlinear operator that represents the Backward Euler nonlinear operator at a
given time step, i.e., residual(t, u_n+1, (u_n+1 - u_n) / dt)
"""
struct BackwardEulerSolverOperator <: NonlinearOperator
  ode_op::ODEOperator
  ode_cache
  tF::Float64
  dt::Float64
  u0::AbstractVector
  u̇F::AbstractVector
end

function Algebra.zero_initial_guess(op::BackwardEulerSolverOperator)
  u = similar(op.u0)
  fill!(u, zero(eltype(u)))
  u
end

function Algebra.allocate_residual(
  op::BackwardEulerSolverOperator,
  u::AbstractVector
)
  allocate_residual(op.ode_op, op.tF, (u, u), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::BackwardEulerSolverOperator,
  u::AbstractVector
)
  uF = u
  u0, u̇F, dt = op.u0, op.u̇F, op.dt
  _discrete_time_derivative!(u̇F, u0, uF, dt)
  residual!(r, op.ode_op, op.tF, (uF, u̇F), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::BackwardEulerSolverOperator,
  u::AbstractVector
)
  allocate_jacobian(op.ode_op, op.tF, (u, u), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::BackwardEulerSolverOperator,
  u::AbstractVector
)
  fillstored!(J, zero(eltype(J)))

  uF = u
  u0, u̇F, dt = op.u0, op.u̇F, op.dt
  _discrete_time_derivative!(u̇F, u0, uF, dt)
  jacobians!(J, op.ode_op, op.tF, (uF, u̇F), (1, 1 / op.dt), op.ode_cache)
end
