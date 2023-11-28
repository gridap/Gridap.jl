"""
    struct ForwardEuler <: ODESolver end

Forward Euler ODE solver.
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
    solver_cache = allocate_solver_cache(ode_op, ode_cache, t0, u0)
  else
    ode_cache, u̇F, nls_cache, solver_cache = cache
  end

  dt = solver.dt
  tF = t0 + dt

  # Update Dirichlet boundary conditions
  ode_cache = update_cache!(ode_cache, ode_op, t0)

  # Create and solve discrete ODE operator
  nl_op = ForwardEulerSolverOperator(
    ode_op, ode_cache, solver_cache,
    t0, dt, u0, u̇F
  )
  nls_cache = solve!(uF, solver.nls, nl_op, nls_cache)

  # Update cache
  cache = (ode_cache, u̇F, nls_cache, solver_cache)

  (uF, tF, cache)
end

function allocate_solver_cache(
  ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  nothing
end

function allocate_solver_cache(
  ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
  r = allocate_residual(ode_op, t, (u, u), ode_cache)
  J, r
end

"""
    ForwardEulerSolverOperator(
      ode_op::ODEOperator, ode_cache, solver_cache,
      t0::Real, dt::Real, u0::AbstractVector, u̇F::AbstractVector
    ) -> Union{ForwardEulerSolverNonlinearOperator, ForwardEulerSolverLinearOperator}

Return the linear or nonlinear Forward Euler operator, depending on the
type of `ODEOperator`.
"""
function ForwardEulerSolverOperator(
  ode_op::ODEOperator, ode_cache, solver_cache,
  t0::Real, dt::Real, u0::AbstractVector, u̇F::AbstractVector
)
  ForwardEulerSolverNonlinearOperator(
    ode_op, ode_cache,
    t0, dt, u0, u̇F
  )
end

function ForwardEulerSolverOperator(
  ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache, solver_cache,
  t0::Real, dt::Real, u0::AbstractVector, u̇F::AbstractVector
)
  ForwardEulerSolverLinearOperator(
    ode_op, ode_cache, solver_cache,
    t0, dt, u0, u̇F
  )
end

"""
    struct ForwardEulerSolverNonlinearOperator <: NonlinearOperator

Forward Euler nonlinear operator at a given time step, i.e.
```math
residual(t_n, u_n, (u_n+1 - u_n) / dt) = 0.
```
"""
struct ForwardEulerSolverNonlinearOperator <: NonlinearOperator
  ode_op::ODEOperator
  ode_cache
  t0::Float64
  dt::Float64
  u0::AbstractVector
  u̇F::AbstractVector
end

function Algebra.zero_initial_guess(op::ForwardEulerSolverNonlinearOperator)
  u = similar(op.u0)
  fill!(u, zero(eltype(u)))
  u
end

function Algebra.allocate_residual(
  op::ForwardEulerSolverNonlinearOperator,
  u::AbstractVector
)
  allocate_residual(op.ode_op, op.t0, (u, u), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ForwardEulerSolverNonlinearOperator,
  u::AbstractVector
)
  u0, u̇F, dt = op.u0, op.u̇F, op.dt
  _discrete_time_derivative!(u̇F, u0, u, dt)
  residual!(r, op.ode_op, op.t0, (u0, u̇F), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ForwardEulerSolverNonlinearOperator,
  u::AbstractVector
)
  allocate_jacobian(op.ode_op, op.t0, (u, u), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ForwardEulerSolverNonlinearOperator,
  u::AbstractVector
)
  fillstored!(J, zero(eltype(J)))

  u0, u̇F, dt = op.u0, op.u̇F, op.dt
  _discrete_time_derivative!(u̇F, u0, u, dt)
  jacobian!(J, op.ode_op, op.t0, (u0, u̇F), 1, inv(dt), op.ode_cache)
end

"""
    struct ForwardEulerSolverLinearOperator <: AffineOperator

Forward Euler linear operator at a given time step, i.e.
`mass(t_n, u_n) * (u_n+1 - u_n) / dt + res(t_n, u_n) = 0`. For simplicity,
we solve for `k_n = (u_n+1 - u_n) / dt` and we recover `u_n+1` from `k_n`,
`u_n` and `dt`.
"""
function ForwardEulerSolverLinearOperator(
  ode_op, ode_cache, solver_cache,
  t0, dt, u0, u̇F
)
  J, r = solver_cache
  AffineOperator(J, r)
end
