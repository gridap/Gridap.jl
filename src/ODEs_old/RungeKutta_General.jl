"""
Abstract Diagonally-Implicit Runge-Kutta ODE solver
"""
abstract type AbstractRungeKutta <: ODESolver end

function num_subsolver(solver::AbstractRungeKutta)
  @abstractmethod
end

function get_subsolver(solver::AbstractRungeKutta, A, i)
  @abstractmethod
end

function get_tableau(solver::AbstractRungeKutta)
  @abstractmethod
end

function allocate_cache(solver::AbstractRungeKutta, ode_op, ode_cache, u0, t0)
  @abstractmethod
end

"""
solve_step!(uf, solver, ode_op, u0, t0, cache)
"""
function solve_step!(
  uf::AbstractVector, solver::AbstractRungeKutta, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real, cache
)
  # Unpack variables
  dt = solver.dt
  tableau = get_tableau(solver)
  A = get_matrix(tableau)
  b = get_weights(tableau)
  c = get_nodes(tableau)
  num_stages = length(b)

  # Allocate cache if not there
  if isnothing(cache)
    ode_cache = allocate_cache(ode_op)
    M = allocate_jacobian(ode_op, t0, u0, ode_cache)
    update_mass_matrix!(M, ode_op, t0, u0, ode_cache)

    solver_cache = allocate_cache(solver, ode_op, ode_cache, u0, t0)
    ks = [similar(u0) for _ in 1:num_stages]

    subsolver_cache = tuple((nothing for _ in 1:num_subsolvers(solver))...)
  else
    ode_cache, M, solver_cache, ks, subsolver_cache = cache
  end

  # Create Runge Kutta ODE operator
  S = typeof(solver)
  rk_op = RungeKuttaOperator{S}(
    ode_op, ode_cache, u0, M,
    solver_cache, 0, t0, dt, ks, A,
  )

  # Solve intermediate stages
  for i in 1:num_stages
    # Update time and caches
    ti = t0 + c[i] * dt
    ode_cache = update_cache!(ode_cache, ode_op, ti)
    update_rk_op!(rk_op, i, ti, ks[i])

    # Solve at stage i
    subsolver, idx = get_subsolver(solver, A, i)
    subcache = subsolver_cache[idx]
    subcache = solve!(uf, subsolver, rk_op, subcache)
    subsolver_cache = Base.setindex(subsolver_cache, subcache, idx)

    # Update stage unknown
    @. ks[i] = uf
  end

  # Update final time
  tf = t0 + dt

  # Take final linear combination
  @. uf = u0
  for i in 1:num_stages
    @. uf = uf + dt * b[i] * ks[i]
  end

  # Update final cache
  cache = (ode_cache, M, solver_cache, ks, subsolver_cache)

  return (uf, tf, cache)
end

abstract type RungeKuttaNonlinearOperator <: NonlinearOperator end

mutable struct RungeKuttaOperator{S} <: RungeKuttaNonlinearOperator
  ode_op::ODEOperator
  ode_cache
  u0::AbstractVector
  M::AbstractMatrix
  solver_cache
  i::Int
  ti::Float64
  dt::Float64
  ks::Vector{AbstractVector}
  A::Matrix
end

function residual!(
  r::AbstractVector, op::RungeKuttaNonlinearOperator,
  u̇::AbstractVector
)
  @abstractmethod
end

function jacobian!(
  J::AbstractMatrix, op::RungeKuttaNonlinearOperator,
  u̇::AbstractVector
)
  @abstractmethod
end

function allocate_residual(op::RungeKuttaNonlinearOperator, u::AbstractVector)
  allocate_residual(op.ode_op, op.ti, u, op.ode_cache)
end

function allocate_jacobian(op::RungeKuttaNonlinearOperator, u::AbstractVector)
  allocate_jacobian(op.ode_op, op.ti, u, op.ode_cache)
end

function zero_initial_guess(op::RungeKuttaNonlinearOperator)
  u̇0 = similar(op.u0)
  fill!(u̇0, zero(eltype(u̇0)))
  u̇0
end

function update_rk_op!(
  op::RungeKuttaNonlinearOperator,
  i::Int, ti::Float64, ki::AbstractVector
)
  op.i = i
  op.ti = ti
  @. op.ks[i] = ki
  op
end

# Redefining solve! function to enforce computation of the jacobian within
# each stage of the Runge-Kutta method when the solver is "LinearSolver".
function solve!(
  x::AbstractVector, ls::LinearSolver,
  op::RungeKuttaNonlinearOperator, cache::Nothing
)
  fill!(x, zero(eltype(x)))
  b = residual(op, x)
  A = jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss, A)
  rmul!(b, -1)
  solve!(x, ns, b)
  LinearSolverCache(A, b, ns)
end

function solve!(
  x::AbstractVector, ls::LinearSolver,
  op::RungeKuttaNonlinearOperator, cache
)
  fill!(x, zero(eltype(x)))
  b = cache.b
  A = cache.A
  ns = cache.ns
  residual!(b, op, x)
  jacobian!(A, op, x)
  numerical_setup!(ns, A)
  rmul!(b, -1)
  solve!(x, ns, b)
  cache
end

# This implementation assumes that the mass matrix is constant w.r.t. u̇
function update_mass_matrix!(
  M::AbstractMatrix, op::ODEOperator,
  t::Float64, u::AbstractVector, cache
)
  z = zero(eltype(M))
  fillstored!(M, z)
  jacobian!(M, op, t, (u, u), 2, 1, cache)
  M
end

"""
Explicit Runge Kutta
"""
struct EXRungeKutta <: AbstractRungeKutta
  ls::LinearSolver
  dt::Float64
  tableau::AbstractButcherTableau{ERK}

  function EXRungeKutta(ls::LinearSolver, dt, name::Symbol)
    tableau = ButcherTableau(name)
    new(ls, dt, tableau)
  end
end

function num_subsolvers(solver::EXRungeKutta)
  1
end

function get_subsolver(solver::EXRungeKutta, A, i)
  solver.ls, 1
end

function get_tableau(solver::EXRungeKutta)
  solver.tableau
end

function allocate_cache(solver::EXRungeKutta, ode_op, ode_cache, u0, t0)
  r_temp = allocate_residual(ode_op, t0, u0, ode_cache)
  u_temp = similar(u0)
  (r_temp, u_temp)
end

function residual!(
  r::AbstractVector, op::RungeKuttaOperator{EXRungeKutta},
  u̇::AbstractVector
)
  ode_op, ode_cache, u0 = op.ode_op, op.ode_cache, op.u0
  i, ti, dt, ks = op.i, op.ti, op.dt, op.ks
  A, solver_cache = op.A, op.solver_cache
  r_temp, u_temp = solver_cache

  # This is independent of u̇ and could be executed in update_rk_op
  u = u_temp
  @. u = u0
  for j = 1:i-1
    @. u = u + dt * A[i, j] * ks[j]
  end

  lhs!(r, ode_op, ti, (u, u̇), ode_cache)
  rhs!(r_temp, ode_op, ti, (u, u̇), ode_cache)
  @. r = r - r_temp

  r
end

function jacobian!(
  J::AbstractMatrix, op::RungeKuttaOperator{EXRungeKutta},
  u̇::AbstractVector
)
  @. J = op.M
  J
end

"""
Diagonally Implicit Runge Kutta
"""
struct RungeKutta <: AbstractRungeKutta
  ls::NonlinearSolver
  nls::NonlinearSolver
  dt::Float64
  tableau::AbstractButcherTableau{DIRK}

  function RungeKutta(ls::NonlinearSolver, nls::NonlinearSolver, dt, name::Symbol)
    tableau = ButcherTableau(name)
    new(ls, nls, dt, tableau)
  end
end

function RungeKutta(nls::NonlinearSolver, dt, name::Symbol)
  RungeKutta(nls, nls, dt, name)
end

function num_subsolvers(solver::RungeKutta)
  2
end

function get_subsolver(solver::RungeKutta, A, i)
  iszero(A[i, i]) ? (solver.ls, 1) : (solver.nls, 2)
end

function get_tableau(solver::RungeKutta)
  solver.tableau
end

function allocate_cache(solver::RungeKutta, ode_op, ode_cache, u0, t0)
  r_temp = allocate_residual(ode_op, t0, u0, ode_cache)
  u_temp = similar(u0)
  (r_temp, u_temp)
end

function residual!(
  r::AbstractVector, op::RungeKuttaOperator{RungeKutta},
  u̇::AbstractVector
)
  ode_op, ode_cache, u0 = op.ode_op, op.ode_cache, op.u0
  i, ti, ks, dt = op.i, op.ti, op.ks, op.dt
  A, solver_cache = op.A, op.solver_cache
  r_temp, u_temp = solver_cache

  u = u_temp
  @. u = u0
  for j in 1:i-1
    @. u = u + dt * A[i, j] * ks[j]
  end
  @. u = u + dt * A[i, i] * u̇

  lhs!(r, ode_op, ti, (u, u̇), ode_cache)
  rhs!(r_temp, ode_op, ti, (u, u̇), ode_cache)
  @. r = r - r_temp

  r
end

function jacobian!(
  J::AbstractMatrix, op::RungeKuttaOperator{RungeKutta},
  u̇::AbstractVector
)
  ode_op, ode_cache, u0 = op.ode_op, op.ode_cache, op.u0
  i, ti, ks, dt = op.i, op.ti, op.ks, op.dt
  A, solver_cache = op.A, op.solver_cache
  _, u_temp = solver_cache

  u = u_temp
  @. u = u0
  for j in 1:i-1
    @. u = u + dt * A[i, j] * ks[j]
  end
  @. u = u + dt * A[i, i] * u̇

  z = zero(eltype(J))
  fillstored!(J, z)

  jacobians!(J, ode_op, ti, (u, u̇), (dt * A[i, i], 1), ode_cache)
  J
end
