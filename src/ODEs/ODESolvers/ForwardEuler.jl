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
    ode_cache = allocate_cache(ode_op, t0, (u0, u0))
    nls_cache = nothing
    solver_cache = allocate_solver_cache(solver, ode_op, ode_cache, t0, u0)
  else
    ode_cache, nls_cache, solver_cache = cache
  end

  dt = solver.dt
  tF = t0 + dt

  # Update Dirichlet boundary conditions
  ode_cache = update_cache!(ode_cache, ode_op, t0)

  # Create and solve discrete ODE operator
  nl_op = ForwardEulerOperator(
    ode_op, ode_cache, solver_cache,
    t0, u0
  )
  v = uF
  nls_cache = solve!(v, solver.nls, nl_op, nls_cache)
  uF = _u_from_v!(uF, u0, dt, v)

  # Update cache
  cache = (ode_cache, nls_cache, solver_cache)

  (uF, tF, cache)
end

function allocate_solver_cache(
  solver::ForwardEuler, ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  nothing
end

function allocate_solver_cache(
  solver::ForwardEuler, ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
  r = allocate_residual(ode_op, t, (u, u), ode_cache)
  (J, r)
end

"""
    ForwardEulerOperator(
      ode_op::ODEOperator, ode_cache, solver_cache,
      t0::Real, u0::AbstractVector
    ) -> Union{ForwardEulerNonlinearOperator, ForwardEulerLinearOperator}

Return the linear or nonlinear Forward Euler operator, depending on the
type of `ODEOperator`.
"""
function ForwardEulerOperator(
  ode_op::ODEOperator, ode_cache, solver_cache,
  t0::Real, u0::AbstractVector
)
  ForwardEulerNonlinearOperator(
    ode_op, ode_cache,
    t0, u0
  )
end

function ForwardEulerOperator(
  ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache, solver_cache,
  t0::Real, u0::AbstractVector
)
  ForwardEulerLinearOperator(
    ode_op, ode_cache, solver_cache,
    t0, u0
  )
end

"""
    struct ForwardEulerNonlinearOperator <: NonlinearOperator

Forward Euler nonlinear operator at a given time step, i.e.
```math
residual(t_n, u_n, (u_n+1 - u_n) / dt) = 0.
```
"""
struct ForwardEulerNonlinearOperator <: NonlinearOperator
  ode_op::ODEOperator
  ode_cache
  t0::Float64
  u0::AbstractVector
end

function Algebra.zero_initial_guess(op::ForwardEulerNonlinearOperator)
  v = similar(op.u0)
  fill!(v, zero(eltype(v)))
  v
end

function Algebra.allocate_residual(
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  allocate_residual(op.ode_op, op.t0, (op.u0, v), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  residual!(r, op.ode_op, op.t0, (op.u0, v), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  allocate_jacobian(op.ode_op, op.t0, (op.u0, v), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ForwardEulerNonlinearOperator,
  v::AbstractVector
)
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, op.ode_op, op.t0, (op.u0, v), 1, 1, op.ode_cache)
end

"""
    ForwardEulerLinearOperator(
      ode_op, ode_cache, solver_cache,
      t0, u0
    ) -> AffineOperator

Forward Euler linear operator at a given time step, i.e.
```math
mass(t_n, u_n) * (u_n+1 - u_n) / dt + res(t_n, u_n) = 0.
```
"""
function ForwardEulerLinearOperator(
  ode_op, ode_cache, solver_cache,
  t0, u0
)
  J, r = solver_cache
  v = u0
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, ode_op, t0, (u0, v), 1, 1, ode_cache)
  residual!(r, ode_op, t0, (u0, v), ode_cache, include_highest=false)
  rmul!(r, -1)

  AffineOperator(J, r)
end
