##############
# RungeKutta #
##############
"""
    abstract type RungeKutta <: ODESolver end

Generic Runge-Kutta ODE solver.

# Mandatory
- [`get_tableau(solver)`](@ref)
- [`get_solver_index(solver, explicit)`](@ref)
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
    get_solver_index(
      solver::RungeKutta, explicit::Bool
    ) -> (NonlinearSolver, Integer)

Depending on whether the stage is explicit or not, return the linear or
nonlinear solver for the discrete ODE operator of the Runge-Kutta ODE solver,
together with its index.
"""
function get_solver_index(solver::RungeKutta, explicit::Bool)
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
  uF, sol_cache = solve_dop!(uF, dop, solver, sol_cache)
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
  J, r = nothing, nothing
  (ulc, J, r)
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
  ulc, J, r = dop_cache
  SequentialRungeKuttaNonlinearOperator(
    ode_op, ode_cache, ulc, vs, tableau, J, r,
    t0, u0, dt, t0, zero(t0)
  )
end

function DiscreteODEOperator(
  solver::ExplicitRungeKutta, ode_op::ODEOperator{<:AbstractMassLinearODE}, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = dop_cache
  SequentialRungeKuttaLinearOperator(
    ode_op, ode_cache, ulc, vs, tableau, J, r,
    t0, u0, dt
  )
end

# RungeKutta interface
function get_tableau(solver::ExplicitRungeKutta)
  solver.tableau
end

function get_solver_index(solver::ExplicitRungeKutta, explicit::Bool)
  (solver.sol, 1)
end

################################
# DiagonallyImplicitRungeKutta #
################################
"""
    struct DiagonallyImplicitRungeKutta <: RungeKutta{DiagonallyImplicitTableau} end

Diagonally-implicit Runge-Kutta ODE solver.
"""
struct DiagonallyImplicitRungeKutta <: RungeKutta{DiagonallyImplicitTableau}
  nlsol::NonlinearSolver
  lsol::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{DiagonallyImplicitTableau}
end

# ODESolver interface
function get_dt(solver::DiagonallyImplicitRungeKutta)
  solver.dt
end

function allocate_dop_cache(
  solver::DiagonallyImplicitRungeKutta,
  ode_op::ODEOperator, ode_cache,
  t::Real, u::AbstractVector
)
  ulc = similar(u)
  J, r = nothing, nothing
  (ulc, J, r)
end

function allocate_dop_cache(
  solver::DiagonallyImplicitRungeKutta,
  ode_op::ODEOperator{MassLinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  ulc = similar(u)
  tableau = get_tableau(solver)
  A = get_matrix(tableau)
  if any(i -> iszero(A[i, i]), axes(A, 2))
    J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
    r = allocate_residual(ode_op, t, (u, u), ode_cache)
  else
    J, r = nothing, nothing
  end
  (ulc, J, r)
end

function allocate_dop_cache(
  solver::DiagonallyImplicitRungeKutta,
  ode_op::ODEOperator{LinearODE}, ode_cache,
  t::Real, u::AbstractVector
)
  ulc = similar(u)
  J = allocate_jacobian(ode_op, t, (u, u), ode_cache)
  r = allocate_residual(ode_op, t, (u, u), ode_cache)
  (ulc, J, r)
end

function allocate_sol_cache(solver::DiagonallyImplicitRungeKutta)
  (nothing, nothing)
end

function DiscreteODEOperator(
  solver::DiagonallyImplicitRungeKutta, ode_op::ODEOperator, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = dop_cache
  SequentialRungeKuttaNonlinearOperator(
    ode_op, ode_cache, ulc, vs, tableau, J, r,
    t0, u0, dt, t0, zero(t0)
  )
end

function DiscreteODEOperator(
  solver::DiagonallyImplicitRungeKutta, ode_op::ODEOperator{LinearODE}, ode_cache,
  dop_cache, t0::Real, u0::AbstractVector, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = dop_cache
  SequentialRungeKuttaLinearOperator(
    ode_op, ode_cache, ulc, vs, tableau, J, r,
    t0, u0, dt
  )
end

# RungeKutta interface
function get_tableau(solver::DiagonallyImplicitRungeKutta)
  solver.tableau
end

function get_solver_index(solver::DiagonallyImplicitRungeKutta, explicit::Bool)
  explicit ? (solver.lsol, 2) : (solver.nlsol, 1)
end

#########################################
# SequentialRungeKuttaNonlinearOperator #
#########################################
"""
    struct SequentialRungeKuttaNonlinearOperator <: DiscreteODEOperator end

Nonlinear operator corresponding to a sequential Runge-Kutta (explicit or
diagonally implicit) scheme:
```math
residual(t_s, u_s, v[s]) = 0,
t_s = t_n + c[s] * dt,
u_s = u_n + ∑_{j < s} A[s, j] * dt * v[j] + A[s, s] * dt * v[j],
``` where `A[s, s]` is zero for an explicit scheme.
"""
mutable struct SequentialRungeKuttaNonlinearOperator <: DiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  ulc::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::Union{Nothing,AbstractMatrix}
  r::Union{Nothing,AbstractVector}
  t0::Real
  u0::AbstractVector
  dt::Real
  t::Real
  Ass::Real
end

function Algebra.allocate_residual(
  op::SequentialRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  !iszero(op.Ass) && axpy!(op.Ass * op.dt, v, u)
  r = allocate_residual(op.ode_op, t, (u, v), op.ode_cache)
  !iszero(op.Ass) && axpy!(-op.Ass * op.dt, v, u)
  r
end

function Algebra.residual!(
  r::AbstractVector,
  op::SequentialRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  !iszero(op.Ass) && axpy!(op.Ass * op.dt, v, u)
  residual!(r, op.ode_op, t, (u, v), op.ode_cache)
  !iszero(op.Ass) && axpy!(-op.Ass * op.dt, v, u)
  r
end

function Algebra.allocate_jacobian(
  op::SequentialRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  !iszero(op.Ass) && axpy!(op.Ass * op.dt, v, u)
  J = allocate_jacobian(op.ode_op, t, (u, v), op.ode_cache)
  !iszero(op.Ass) && axpy!(-op.Ass * op.dt, v, u)
  J
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::SequentialRungeKuttaNonlinearOperator,
  v::AbstractVector
)
  t, u = op.t, op.ulc
  !iszero(op.Ass) && axpy!(op.Ass * op.dt, v, u)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, op.ode_op, t, (u, v), (op.Ass * op.dt, 1), op.ode_cache)
  !iszero(op.Ass) && axpy!(-op.Ass * op.dt, v, u)
  J
end

function solve_dop!(
  uF::AbstractVector,
  op::SequentialRungeKuttaNonlinearOperator,
  solver::RungeKutta, caches
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  ulc, vs, tableau = op.ulc, op.vs, op.tableau
  t0, u0, dt, t = op.t0, op.u0, op.dt, op.t

  ismasslinear = ODEOperatorType(ode_op) <: AbstractMassLinearODE
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

    # Update operator state
    op.t = ts
    op.Ass = A[s, s]

    # Solve stage
    explicit = iszero(A[s, s])
    sol, isol = get_solver_index(solver, explicit)
    cache = caches[isol]

    v = vs[s]
    fill!(v, zero(eltype(v)))

    if explicit && ismasslinear
      J, r = op.J, op.r

      fillstored!(J, zero(eltype(J)))
      jacobians!(J, ode_op, ts, (u, v), (A[s, s] * dt, 1), ode_cache)
      residual!(r, ode_op, ts, (u, v), ode_cache, include_highest=false)
      rmul!(r, -1)

      _op = SequentialRungeKuttaLinearOperator(
        ode_op, ode_cache, ulc, vs, tableau, J, r,
        t0, u0, dt
      )
    else
      _op = op
    end

    cache = solve!(v, sol, _op, cache)
    caches = Base.setindex(caches, cache, isol)
  end

  # Take final linear combination
  copy!(uF, u0)
  for s in 1:num_stages
    coef = b[s]
    if !iszero(coef)
      axpy!(coef * dt, vs[s], uF)
    end
  end

  (uF, caches)
end

######################################
# SequentialRungeKuttaLinearOperator #
######################################
"""
    struct SequentialRungeKuttaLinearOperator <: LinearDiscreteODEOperator end

Linear operator corresponding to a sequential Runge-Kutta (explicit or
diagonally implicit) scheme:
```math
residual(t_s, u_s, v[s]) = mass(t_s, u_s) v[s] + res(t_s, u_s) = 0,
t_s = t_n + c[s] * dt,
u_s = u_n + ∑_{j < s} A[s, j] * dt * v[j] + A[s, s] * dt * v[s],
``` where `A[s, s]` is zero for an explicit scheme.
"""
struct SequentialRungeKuttaLinearOperator <: LinearDiscreteODEOperator
  ode_op::ODEOperator
  ode_cache
  ulc::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  u0::AbstractVector
  dt::Real
end

Algebra.get_matrix(op::SequentialRungeKuttaLinearOperator) = op.J

Algebra.get_vector(op::SequentialRungeKuttaLinearOperator) = op.r

function solve_dop!(
  uF::AbstractVector,
  op::SequentialRungeKuttaLinearOperator,
  solver::RungeKutta, caches
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  ulc, vs, tableau = op.ulc, op.vs, op.tableau
  J, r = op.J, op.r
  t0, u0, dt = op.t0, op.u0, op.dt

  explicit = true
  sol, isol = get_solver_index(solver, explicit)
  cache = caches[isol]

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
    jacobians!(J, ode_op, ts, (u, v), (A[s, s] * dt, 1), ode_cache)
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

  caches = Base.setindex(caches, cache, isol)

  (uF, caches)
end
