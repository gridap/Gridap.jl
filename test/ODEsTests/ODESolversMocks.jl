using LinearAlgebra
using LinearAlgebra: fillstored!

using Gridap
using Gridap.Algebra
using Gridap.ODEs
using Gridap.ODEs: _u_from_v!

################
# NLSolverMock #
################
struct NLSolverMock <: NonlinearSolver
end

function Algebra.solve!(
  x::AbstractVector, nls::NLSolverMock,
  op::NonlinearOperator, cache::Nothing
)
  r = residual(op, x)
  J = jacobian(op, x)

  dx = -J \ r
  axpy!(1, dx, x)

  cache = (r, J, dx)
  cache
end

function Algebra.solve!(
  x::AbstractVector, nls::NLSolverMock,
  op::NonlinearOperator, cache
)
  r, J, dx = cache

  residual!(r, op, x)
  jacobian!(J, op, x)

  dx = -J \ r
  axpy!(1, dx, x)

  cache
end

#################
# ODESolverMock #
#################
struct ODESolverMock <: ODESolver
  sol::NLSolverMock
  dt::Float64
end

function ODEs.get_dt(solver::ODESolverMock)
  solver.dt
end

function ODEs.allocate_dop_cache(
  solver::ODESolverMock,
  ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  nothing
end

function ODEs.allocate_sol_cache(solver::ODESolverMock)
  nothing
end

function ODEs.DiscreteODEOperator(
  solver::ODESolverMock, ode_op::ODEOperator, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real
)
  MockDiscreteODEOperator(ode_op, ode_cache, t0, u0, dt)
end

function ODEs.solve_step!(
  uF::AbstractVector,
  solver::ODESolverMock, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real, cache
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
  sol_cache = solve_dop!(uF, dop, solver.sol, sol_cache)
  tF = t0 + dt

  # Update cache
  cache = (ode_cache, dop_cache, sol_cache)

  (uF, tF, cache)
end

###########################
# MockDiscreteODEOperator #
###########################
struct MockDiscreteODEOperator <: DiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  t0::Real
  u0::AbstractVector
  dt::Real
end

function Algebra.allocate_residual(
  op::MockDiscreteODEOperator,
  v::AbstractVector
)
  t, u0 = op.t0 + op.dt, op.u0
  u = u0 .+ dt .* v
  allocate_residual(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector, op::MockDiscreteODEOperator, v::AbstractVector
)
  t, u0 = op.t0 + op.dt, op.u0
  u = u0 .+ dt .* v
  residual!(r, op.ode_op, t, (u, v), op.ode_cache)
  r
end

function Algebra.allocate_jacobian(
  op::MockDiscreteODEOperator, v::AbstractVector
)
  t, u0 = op.t0 + op.dt, op.u0
  u = u0 .+ dt .* v
  allocate_jacobian(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::MockDiscreteODEOperator, v::AbstractVector
)
  t, u0 = op.t0 + op.dt, op.u0
  u = u0 .+ dt .* v
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, op.ode_op, t, (u, v), (op.dt, 1), op.ode_cache)
  J
end

function ODEs.solve_dop!(
  uF::AbstractVector,
  op::MockDiscreteODEOperator, sol::NonlinearSolver, cache
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  t0, u0, dt = op.t0, op.u0, op.dt

  update_cache!(ode_cache, ode_op, t0 + dt)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  cache = solve!(vF, sol, op, cache)
  _u_from_v!(uF, u0, dt, vF)

  cache
end
