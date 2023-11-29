##############
# RungeKutta #
##############
"""
    abstract type RungeKutta <: ODESolver end

Generic Runge-Kutta ODE solver.

# Mandatory
- [`get_tableau(solver)`](@ref)
- [`get_sols(solver)`](@ref)
"""
abstract type RungeKutta{T<:TableauType} <: ODESolver end

function RungeKutta(
  nls::NonlinearSolver, ls::NonlinearSolver,
  dt::Real, name::Symbol
)
  tableau = ButcherTableau(name)
  type = TableauType(tableau)
  if type == ExplicitTableau
    ExplicitRungeKutta(nls, dt, tableau)
  elseif type == DiagonallyImplicitTableau
    DiagonallyImplicitRungeKutta(nls, ls, dt, tableau)
    # elseif type == FullyImplicitTableau
    #   FullyImplicitRungeKutta(nls, ls, dt, tableau)
  end
end

function RungeKutta(ls::NonlinearSolver, dt::Real, name::Symbol)
  RungeKutta(ls, ls, dt, name)
end

"""
    get_tableau(solver::RungeKutta) -> Integer

Return the Butcher tableau of the Runge-Kutta ODE solver.
"""
function get_tableau(solver::RungeKutta)
  @abstractmethod
end

"""
    get_sols(solver::RungeKutta) -> Integer

Return the solvers for the discrete ODE operator of the Runge-Kutta ODE solver.
"""
function get_sols(solver::RungeKutta)
  @abstractmethod
end

function solve_step!(
  uF::AbstractVector,
  solver::RungeKutta, ode_op::ODEOperator,
  u0::AbstractVector, t0::Real,
  cache
)
  # Unpack solver
  dt = get_dt(solver)
  tableau = get_tableau(solver)
  num_stages = length(get_nodes(tableau))

  # Allocate or unpack cache
  if isnothing(cache)
    ode_cache = allocate_cache(ode_op, t0, (u0, u0))
    dop_cache = allocate_dop_cache(solver, ode_op, ode_cache, t0, u0)
    sol_cache = allocate_sol_cache(solver)
    vs = [similar(u0) for _ in 1:num_stages]
  else
    ode_cache, dop_cache, sol_cache, vs = cache
  end

  # Create discrete ODE operator
  dop = DiscreteODEOperator(
    solver, ode_op, ode_cache, dop_cache,
    t0, u0, dt, vs, tableau
  )

  # Solve discrete ODE operator
  sol_cache = solve_dop!(uF, dop, get_sols(solver), sol_cache)
  tF = t0 + dt

  # Update cache
  cache = (ode_cache, dop_cache, sol_cache, vs)

  (uF, tF, cache)
end

######################
# ExplicitRungeKutta #
######################
"""
    struct ExplicitRungeKutta <: RungeKutta{ExplicitTableau} end

Explicit Runge-Kutta ODE solver.
"""
struct ExplicitRungeKutta <: RungeKutta{ExplicitTableau}
  sol::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{ExplicitTableau}

  function ExplicitRungeKutta(
    sol::NonlinearSolver, dt::Real,
    tableau::AbstractTableau{ExplicitTableau}
  )
    new(sol, dt, tableau)
  end
end

# ODESolver interface
function get_dt(solver::ExplicitRungeKutta)
  solver.dt
end

function allocate_dop_cache(
  solver::ExplicitRungeKutta,
  ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  ulc = similar(u)
  (ulc,)
end

function allocate_dop_cache(
  solver::ExplicitRungeKutta,
  ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  ulc = similar(u)
  J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
  r = allocate_residual(ode_op, t, (u, u), ode_cache)
  (ulc, J, r)
end

function allocate_sol_cache(solver::ExplicitRungeKutta)
  (nothing,)
end

function DiscreteODEOperator(
  solver::ExplicitRungeKutta, ode_op::ODEOperator, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, = dop_cache
  ExplicitRungeKuttaNonlinearOperator(
    ode_op, ode_cache, ulc, vs, tableau,
    t0, u0, dt, t0
  )
end

function DiscreteODEOperator(
  solver::ExplicitRungeKutta, ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = dop_cache
  ExplicitRungeKuttaLinearOperator(
    ode_op, ode_cache, ulc, vs, tableau, J, r,
    t0, u0, dt
  )
end

# RungeKutta interface
function get_tableau(solver::ExplicitRungeKutta)
  solver.tableau
end

function get_sols(solver::ExplicitRungeKutta)
  (solver.sol,)
end

#######################################
# ExplicitRungeKuttaNonlinearOperator #
#######################################
"""
    struct ExplicitRungeKuttaNonlinearOperator <: DiscreteODEOperator end

Nonlinear operator corresponding to an explicit Runge-Kutta scheme:
```math
residual(t_s, u_s, v[s]) = 0,
t_s = t_n + c[s] * dt,
u_s = u_n + ∑_{j < s} A[s, j] * dt * v[j]
```
"""
mutable struct ExplicitRungeKuttaNonlinearOperator <: DiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  ulc::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau{ExplicitTableau}
  t0::Real
  u0::AbstractVector
  dt::Real
  t::Real
end

function Algebra.allocate_residual(
  op::ExplicitRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  allocate_residual(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::ExplicitRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  residual!(r, op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::ExplicitRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  allocate_jacobian(op.ode_op, t, (u, v), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::ExplicitRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, op.ode_op, t, (u, v), 1, 1, op.ode_cache)
end

function solve_dop!(
  uF::AbstractVector,
  op::ExplicitRungeKuttaNonlinearOperator,
  sols::Tuple{Vararg{NonlinearSolver}}, caches
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  ulc, vs, tableau = op.ulc, op.vs, op.tableau
  t0, u0, dt, t = op.t0, op.u0, op.dt, op.t
  sol, cache = first(sols), first(caches)

  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  num_stages = length(c)

  # Solve stages
  for s in 1:num_stages
    ts = t0 + c[s] * dt
    update_cache!(ode_cache, ode_op, ts)
    op.t = ts

    # Take linear combination of previous stages
    u = ulc
    copy!(u, u0)
    for j in 1:s-1
      coef = A[s, j]
      if !iszero(coef)
        axpy!(coef * dt, vs[j], u)
      end
    end

    # Solve stage
    v = vs[s]
    fill!(v, zero(eltype(v)))
    cache = solve!(v, sol, op, cache)
  end

  # Take final linear combination
  copy!(uF, u0)
  for s in 1:num_stages
    coef = b[s]
    if !iszero(coef)
      axpy!(coef * dt, vs[s], uF)
    end
  end

  caches = (cache,)
  caches
end

####################################
# ExplicitRungeKuttaLinearOperator #
####################################
"""
    struct ExplicitRungeKuttaLinearOperator <: LinearDiscreteODEOperator end

Linear operator corresponding to an explicit Runge-Kutta scheme:
```math
residual(t_s, u_s, v[s]) = mass(t_s, u_s) v[s] + res(t_s, u_s) = 0,
t_s = t_n + c[s] * dt,
u_s = u_n + ∑_{j < s} A[s, j] * dt * v[j]
```
"""
struct ExplicitRungeKuttaLinearOperator <: LinearDiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  ulc::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau{ExplicitTableau}
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  u0::AbstractVector
  dt::Real
end

Algebra.get_matrix(op::ExplicitRungeKuttaLinearOperator) = op.J

Algebra.get_vector(op::ExplicitRungeKuttaLinearOperator) = op.r

function solve_dop!(
  uF::AbstractVector,
  op::ExplicitRungeKuttaLinearOperator,
  sols::Tuple{Vararg{NonlinearSolver}}, caches
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  ulc, vs, tableau = op.ulc, op.vs, op.tableau
  J, r = op.J, op.r
  t0, u0, dt = op.t0, op.u0, op.dt
  sol, cache = first(sols), first(caches)

  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  num_stages = length(c)

  # Solve stages
  for s in 1:num_stages
    ts = t0 + c[s] * dt
    update_cache!(ode_cache, ode_op, ts)

    # Take linear combination of previous stages
    u = ulc
    copy!(u, u0)
    for j in 1:s-1
      coef = A[s, j]
      if !iszero(coef)
        axpy!(coef * dt, vs[j], u)
      end
    end
    v = vs[s]

    # Update jacobian and residual
    fillstored!(J, zero(eltype(J)))
    jacobian!(J, ode_op, ts, (u, v), 1, 1, ode_cache)
    residual!(r, ode_op, ts, (u, v), ode_cache, include_highest=false)
    rmul!(r, -1)

    # Solve stage
    fill!(v, zero(eltype(v)))
    cache = solve!(v, sol, op, cache)
  end

  # Take final linear combination
  copy!(uF, u0)
  for s in 1:num_stages
    coef = b[s]
    if !iszero(coef)
      axpy!(coef * dt, vs[s], uF)
    end
  end

  caches = (cache,)
  caches
end
