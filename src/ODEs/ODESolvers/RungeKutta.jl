##############
# RungeKutta #
##############
"""
    abstract type RungeKutta <: ODESolver end

Generic Runge-Kutta ODE solver.

# Mandatory
- [`get_tableau(odeslvr)`](@ref)
- [`get_solver_index(odeslvr, explicit)`](@ref)
"""
abstract type RungeKutta{T<:TableauType} <: ODESolver end

function RungeKutta(
  disop_nl::NonlinearSolver, disop_l::NonlinearSolver,
  dt::Real, name::Symbol
)
  tableau = ButcherTableau(name)
  type = TableauType(tableau)
  if type == ExplicitTableau
    ExplicitRungeKutta(disop_nl, dt, tableau)
  elseif type == DiagonallyImplicitTableau
    DiagonallyImplicitRungeKutta(disop_nl, disop_l, dt, tableau)
    # elseif type == FullyImplicitTableau
    #   FullyImplicitRungeKutta(disop_nl, disop_l, dt, tableau)
  end
end

function RungeKutta(disop_nl::NonlinearSolver, dt::Real, name::Symbol)
  RungeKutta(disop_nl, disop_nl, dt, name)
end

"""
    get_tableau(odeslvr::RungeKutta) -> Integer

Return the Butcher tableau of the Runge-Kutta ODE solver.
"""
function get_tableau(odeslvr::RungeKutta)
  @abstractmethod
end

"""
    get_solver_index(
      odeslvr::RungeKutta, explicit::Bool
    ) -> (NonlinearSolver, Integer)

Depending on whether the stage is explicit or not, return the linear or
nonlinear solver for the discrete ODE operator of the Runge-Kutta ODE solver,
together with its index.
"""
function get_solver_index(odeslvr::RungeKutta, explicit::Bool)
  @abstractmethod
end

function solve_step!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, odeop::ODEOperator,
  us0::NTuple{1,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, = us0
  dt = get_dt(odeslvr)
  tableau = get_tableau(odeslvr)
  num_stages = length(get_nodes(tableau))

  # Allocate or unpack cache
  if isnothing(cache)
    odeopcache = allocate_odeopcache(odeop, t0, (u0, u0))
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, u0)
    disslvrcache = allocate_disslvrcache(odeslvr)
    vs = [similar(u0) for _ in 1:num_stages]
  else
    odeopcache, disopcache, disslvrcache, vs = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop,
    odeopcache, disopcache,
    t0, us0, dt,
    vs, tableau
  )

  # Solve discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache, vs)

  (usF, tF, cache)
end

######################
# ExplicitRungeKutta #
######################
"""
    struct ExplicitRungeKutta <: RungeKutta{ExplicitTableau} end

Explicit Runge-Kutta ODE solver.
"""
struct ExplicitRungeKutta <: RungeKutta{ExplicitTableau}
  disslvr::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{ExplicitTableau}
end

# ODESolver interface
function get_dt(odeslvr::ExplicitRungeKutta)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::ExplicitRungeKutta,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  ulc = similar(x)
  J, r = nothing, nothing
  (ulc, J, r)
end

function allocate_disopcache(
  odeslvr::ExplicitRungeKutta,
  odeop::ODEOperator{<:AbstractMassLinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  ulc = similar(x)
  J = allocate_jacobian(odeop, t, (x, x), odeopcache)
  r = allocate_residual(odeop, t, (x, x), odeopcache)
  (ulc, J, r)
end

function allocate_disslvrcache(odeslvr::ExplicitRungeKutta)
  (nothing,)
end

function DiscreteODEOperator(
  odeslvr::ExplicitRungeKutta, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = disopcache
  SequentialRungeKuttaNonlinearOperator(
    odeop, odeopcache, ulc, vs, tableau, J, r,
    t0, us0, dt, t0, zero(t0)
  )
end

function DiscreteODEOperator(
  odeslvr::ExplicitRungeKutta, odeop::ODEOperator{<:AbstractMassLinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = disopcache
  SequentialRungeKuttaLinearOperator(
    odeop, odeopcache, ulc, vs, tableau, J, r,
    t0, us0, dt
  )
end

# RungeKutta interface
function get_tableau(odeslvr::ExplicitRungeKutta)
  odeslvr.tableau
end

function get_solver_index(odeslvr::ExplicitRungeKutta, explicit::Bool)
  (odeslvr.disslvr, 1)
end

################################
# DiagonallyImplicitRungeKutta #
################################
"""
    struct DiagonallyImplicitRungeKutta <: RungeKutta{DiagonallyImplicitTableau} end

Diagonally-implicit Runge-Kutta ODE solver.
"""
struct DiagonallyImplicitRungeKutta <: RungeKutta{DiagonallyImplicitTableau}
  disslvr_nl::NonlinearSolver
  disslvr_l::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{DiagonallyImplicitTableau}
end

# ODESolver interface
function get_dt(odeslvr::DiagonallyImplicitRungeKutta)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::DiagonallyImplicitRungeKutta,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  ulc = similar(x)
  J, r = nothing, nothing
  (ulc, J, r)
end

function allocate_disopcache(
  odeslvr::DiagonallyImplicitRungeKutta,
  odeop::ODEOperator{MassLinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  ulc = similar(x)
  tableau = get_tableau(odeslvr)
  A = get_matrix(tableau)
  if any(i -> iszero(A[i, i]), axes(A, 2))
    J = allocate_jacobian(odeop, t, (x, x), odeopcache)
    r = allocate_residual(odeop, t, (x, x), odeopcache)
  else
    J, r = nothing, nothing
  end
  (ulc, J, r)
end

function allocate_disopcache(
  odeslvr::DiagonallyImplicitRungeKutta,
  odeop::ODEOperator{LinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  ulc = similar(x)
  J = allocate_jacobian(odeop, t, (x, x), odeopcache)
  r = allocate_residual(odeop, t, (x, x), odeopcache)
  (ulc, J, r)
end

function allocate_disslvrcache(odeslvr::DiagonallyImplicitRungeKutta)
  (nothing, nothing)
end

function DiscreteODEOperator(
  odeslvr::DiagonallyImplicitRungeKutta, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = disopcache
  SequentialRungeKuttaNonlinearOperator(
    odeop, odeopcache, ulc, vs, tableau, J, r,
    t0, us0, dt, t0, zero(t0)
  )
end

function DiscreteODEOperator(
  odeslvr::DiagonallyImplicitRungeKutta, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ulc, J, r = disopcache
  SequentialRungeKuttaLinearOperator(
    odeop, odeopcache, ulc, vs, tableau, J, r,
    t0, us0, dt
  )
end

# RungeKutta interface
function get_tableau(odeslvr::DiagonallyImplicitRungeKutta)
  odeslvr.tableau
end

function get_solver_index(odeslvr::DiagonallyImplicitRungeKutta, explicit::Bool)
  explicit ? (odeslvr.disslvr_l, 2) : (odeslvr.disslvr_nl, 1)
end

#########################################
# SequentialRungeKuttaNonlinearOperator #
#########################################
"""
    struct SequentialRungeKuttaNonlinearOperator <: DiscreteODEOperator end

Nonlinear operator corresponding to a sequential Runge-Kutta (explicit or
diagonally implicit) scheme:
```math
residual(t_s, u_s, x[s]) = 0,
t_s = t_n + c[s] * dt,
u_s = u_n + ∑_{j < s} A[s, j] * dt * x[j] + A[s, s] * dt * x[j],
u_n+1 = u_n + ∑_{i} b_i * dt * x[i]
``` where `A[s, s]` is zero for an explicit scheme.
"""
mutable struct SequentialRungeKuttaNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  ulc::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::Union{Nothing,AbstractMatrix}
  r::Union{Nothing,AbstractVector}
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  t::Real
  Ass::Real
end

function Algebra.allocate_residual(
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  t, dt, u = disop.t, disop.dt, disop.ulc
  !iszero(disop.Ass) && axpy!(disop.Ass * dt, x, u)
  r = allocate_residual(disop.odeop, t, (u, x), disop.odeopcache)
  !iszero(disop.Ass) && axpy!(-disop.Ass * dt, x, u)
  r
end

function Algebra.residual!(
  r::AbstractVector,
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  t, dt, u = disop.t, disop.dt, disop.ulc
  !iszero(disop.Ass) && axpy!(disop.Ass * dt, x, u)
  residual!(r, disop.odeop, t, (u, x), disop.odeopcache)
  !iszero(disop.Ass) && axpy!(-disop.Ass * dt, x, u)
  r
end

function Algebra.allocate_jacobian(
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  t, dt, u = disop.t, disop.dt, disop.ulc
  !iszero(disop.Ass) && axpy!(disop.Ass * dt, x, u)
  J = allocate_jacobian(disop.odeop, t, (u, x), disop.odeopcache)
  !iszero(disop.Ass) && axpy!(-disop.Ass * dt, x, u)
  J
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  t, dt, u = disop.t, disop.dt, disop.ulc
  !iszero(disop.Ass) && axpy!(disop.Ass * dt, x, u)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, t, (u, x), (disop.Ass * dt, 1), disop.odeopcache)
  !iszero(disop.Ass) && axpy!(-disop.Ass * dt, x, u)
  J
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, disop::SequentialRungeKuttaNonlinearOperator,
  disslvrcaches
)
  uF, = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  ulc, vs, tableau = disop.ulc, disop.vs, disop.tableau
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0

  ismasslinear = ODEOperatorType(odeop) <: AbstractMassLinearODE
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  num_stages = length(c)

  # Solve stages
  for s in 1:num_stages
    ts = t0 + c[s] * dt
    update_odeopcache!(odeopcache, odeop, ts)

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
    disop.t = ts
    disop.Ass = A[s, s]

    # Solve stage
    explicit = iszero(A[s, s])
    disslvr, islvr = get_solver_index(odeslvr, explicit)
    disslvrcache = disslvrcaches[islvr]

    x = vs[s]
    fill!(x, zero(eltype(x)))

    if explicit && ismasslinear
      J, r = disop.J, disop.r

      fillstored!(J, zero(eltype(J)))
      jacobians!(J, odeop, ts, (u, x), (A[s, s] * dt, 1), odeopcache)
      residual!(r, odeop, ts, (u, x), odeopcache, include_highest=false)
      rmul!(r, -1)

      _op = SequentialRungeKuttaLinearOperator(
        odeop, odeopcache, ulc, vs, tableau, J, r,
        t0, disop.us0, dt
      )
    else
      _op = disop
    end

    disslvrcache = solve!(x, disslvr, _op, disslvrcache)
    disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)
  end

  # Take final linear combination
  copy!(uF, u0)
  for s in 1:num_stages
    coef = b[s]
    if !iszero(coef)
      axpy!(coef * dt, vs[s], uF)
    end
  end

  usF = (uF,)
  (usF, disslvrcaches)
end

######################################
# SequentialRungeKuttaLinearOperator #
######################################
"""
    struct SequentialRungeKuttaLinearOperator <: LinearDiscreteODEOperator end

Linear operator corresponding to a sequential Runge-Kutta (explicit or
diagonally implicit) scheme:
```math
residual(t_s, u_s, x[s]) = mass(t_s, u_s) x[s] + res(t_s, u_s) = 0,
t_s = t_n + c[s] * dt,
u_s = u_n + ∑_{j < s} A[s, j] * dt * x[j] + A[s, s] * dt * x[s],
``` where `A[s, s]` is zero for an explicit scheme.
"""
struct SequentialRungeKuttaLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  ulc::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
end

Algebra.get_matrix(disop::SequentialRungeKuttaLinearOperator) = disop.J

Algebra.get_vector(disop::SequentialRungeKuttaLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, disop::SequentialRungeKuttaLinearOperator,
  disslvrcaches
)
  uF, = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  ulc, vs, tableau = disop.ulc, disop.vs, disop.tableau
  J, r = disop.J, disop.r
  t0, dt, (u0,) = disop.t0, disop.dt, disop.us0

  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  num_stages = length(c)

  # Solve stages
  for s in 1:num_stages
    ts = t0 + c[s] * dt
    update_odeopcache!(odeopcache, odeop, ts)

    # Take linear combination of previous stages
    u = ulc
    copy!(u, u0)
    for j in 1:s-1
      coef = A[s, j]
      if !iszero(coef)
        axpy!(coef * dt, vs[j], u)
      end
    end
    x = vs[s]

    # Update jacobian and residual
    fillstored!(J, zero(eltype(J)))
    jacobians!(J, odeop, ts, (u, x), (A[s, s] * dt, 1), odeopcache)
    residual!(r, odeop, ts, (u, x), odeopcache, include_highest=false)
    rmul!(r, -1)

    # Solve stage
    fill!(x, zero(eltype(x)))
    disslvrcache = solve!(x, disslvr, disop, disslvrcache)
  end

  # Take final linear combination
  copy!(uF, u0)
  for s in 1:num_stages
    coef = b[s]
    if !iszero(coef)
      axpy!(coef * dt, vs[s], uF)
    end
  end

  usF = (uF,)
  disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)
  (usF, disslvrcaches)
end
