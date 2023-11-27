using LinearAlgebra

using Gridap
using Gridap.Algebra
using Gridap.ODEs

################
# NLSolverMock #
################
struct NLSolverMock <: NonlinearSolver
end

function Gridap.Algebra.solve!(
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

function Gridap.Algebra.solve!(
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
  nls::NLSolverMock
  dt::Float64
end

function Gridap.ODEs.solve_step!(
  uF::AbstractVector,
  solver::ODESolverMock, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real, cache
)
  dt = solver.dt
  tF = t0 + dt

  # Create or retrieve cache
  if isnothing(cache)
    ode_cache = allocate_cache(ode_op)
    nl_cache = nothing
  else
    ode_cache, nl_cache = cache
  end

  # Update ODE cache
  update_cache!(ode_cache, ode_op, tF)

  # Create and solve ODE operator
  nl_op = OperatorMock(ode_op, u0, dt, tF, ode_cache)
  nl_cache = solve!(uF, solver.nls, nl_op, nl_cache)

  # Update cache
  cache = (ode_cache, nl_cache)

  (uF, tF, cache)
end

################
# OperatorMock #
################
struct OperatorMock <: NonlinearOperator
  ode_op::ODEOperator
  u0::AbstractVector
  dt::Float64
  tF::Float64
  ode_cache
end

function Gridap.Algebra.allocate_residual(op::OperatorMock, u::AbstractVector)
  allocate_residual(op.ode_op, op.tF, (u, u), op.ode_cache)
end

function Gridap.Algebra.residual!(
  r::AbstractVector, op::OperatorMock, u::AbstractVector
)
  u̇ = (u - op.u0) ./ op.dt
  residual!(r, op.ode_op, op.tF, (u, u̇), op.ode_cache)
  r
end

function Gridap.Algebra.allocate_jacobian(op::OperatorMock, u::AbstractVector)
  allocate_jacobian(op.ode_op, op.tF, (u, u), op.ode_cache)
end

function Gridap.Algebra.jacobian!(
  J::AbstractMatrix, op::OperatorMock, u::AbstractVector
)
  u̇ = (u - op.u0) ./ op.dt
  fill!(J, 0)
  jacobians!(J, op.ode_op, op.tF, (u, u̇), (1, inv(op.dt)), op.ode_cache)
  J
end

function Gridap.Algebra.zero_initial_guess(op::OperatorMock)
  u = similar(op.u0)
  fill!(u, zero(eltype(u)))
  u
end
