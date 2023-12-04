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
    EXRungeKutta(disop_nl, dt, tableau)
  elseif type == DiagonallyImplicitTableau
    DIMRungeKutta(disop_nl, disop_l, dt, tableau)
  elseif type == ImplicitExplicitTableau
    IMEXRungeKutta(disop_nl, disop_l, dt, tableau)
    # elseif type == FullyImplicitTableau
    #   FIMRungeKutta(disop_nl, disop_l, dt, tableau)
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
    us = (u0, u0)
    odeopcache = allocate_odeopcache(odeop, t0, us)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, u0)
    disslvrcache = allocate_disslvrcache(odeslvr)
    vs = [zero(u0) for _ in 1:num_stages]
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

  # Solve the discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache, vs)

  (usF, tF, cache)
end

################
# EXRungeKutta #
################
"""
    struct EXRungeKutta <: RungeKutta{ExplicitTableau} end

Explicit Runge-Kutta ODE solver.
"""
struct EXRungeKutta <: RungeKutta{ExplicitTableau}
  disslvr::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{ExplicitTableau}
end

# ODESolver interface
function get_dt(odeslvr::EXRungeKutta)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::EXRungeKutta,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  ui = zero(x)
  J, r = nothing, nothing
  (ui, J, r)
end

function allocate_disopcache(
  odeslvr::EXRungeKutta,
  odeop::ODEOperator{<:AbstractMassLinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  ui = zero(x)
  us = (x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (ui, J, r)
end

function allocate_disslvrcache(odeslvr::EXRungeKutta)
  (nothing,)
end

function DiscreteODEOperator(
  odeslvr::EXRungeKutta, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  ti, aii = t0, zero(t0)
  SequentialRungeKuttaNonlinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    ui, vs, tableau, ti, aii, J, r
  )
end

function DiscreteODEOperator(
  odeslvr::EXRungeKutta, odeop::ODEOperator{<:AbstractMassLinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  SequentialRungeKuttaLinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    ui, vs, tableau, J, r,
  )
end

# RungeKutta interface
function get_tableau(odeslvr::EXRungeKutta)
  odeslvr.tableau
end

function get_solver_index(odeslvr::EXRungeKutta, explicit::Bool)
  (odeslvr.disslvr, 1)
end

#################
# DIMRungeKutta #
#################
"""
    struct DIMRungeKutta <: RungeKutta{DiagonallyImplicitTableau} end

Diagonally-implicit Runge-Kutta ODE solver.
"""
struct DIMRungeKutta <: RungeKutta{DiagonallyImplicitTableau}
  disslvr_nl::NonlinearSolver
  disslvr_l::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{DiagonallyImplicitTableau}
end

# ODESolver interface
function get_dt(odeslvr::DIMRungeKutta)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::DIMRungeKutta,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  ui = zero(x)
  J, r = nothing, nothing
  (ui, J, r)
end

function allocate_disopcache(
  odeslvr::DIMRungeKutta,
  odeop::ODEOperator{MassLinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  ui = zero(x)
  tableau = get_tableau(odeslvr)
  A = get_matrix(tableau)
  if any(i -> iszero(A[i, i]), axes(A, 2))
    us = (x, x)
    J = allocate_jacobian(odeop, t, us, odeopcache)
    r = allocate_residual(odeop, t, us, odeopcache)
  else
    J, r = nothing, nothing
  end
  (ui, J, r)
end

function allocate_disopcache(
  odeslvr::DIMRungeKutta,
  odeop::ODEOperator{LinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  ui = zero(x)
  us = (x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (ui, J, r)
end

function allocate_disslvrcache(odeslvr::DIMRungeKutta)
  (nothing, nothing)
end

function DiscreteODEOperator(
  odeslvr::DIMRungeKutta, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  ti, aii = t0, zero(t0)
  SequentialRungeKuttaNonlinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    ui, vs, tableau, ti, aii, J, r
  )
end

function DiscreteODEOperator(
  odeslvr::DIMRungeKutta, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  vs::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  SequentialRungeKuttaLinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    ui, vs, tableau, J, r
  )
end

# RungeKutta interface
function get_tableau(odeslvr::DIMRungeKutta)
  odeslvr.tableau
end

function get_solver_index(odeslvr::DIMRungeKutta, explicit::Bool)
  explicit ? (odeslvr.disslvr_l, 2) : (odeslvr.disslvr_nl, 1)
end

######################
# Nonlinear operator #
######################
"""
    struct SequentialRungeKuttaNonlinearOperator <: DiscreteODEOperator end

Nonlinear operator corresponding to a sequential Runge-Kutta (explicit or
diagonally implicit) scheme:
```math
residual(ti, ui, vi) = 0,

ti = t_n + c[i] * dt
ui = u_n + dt * ∑_{1 ≤ j < i} A[i, j] * vj + dt * A[i, i] * dt * x
vi = x,

u_(n+1) = u_n + dt * ∑_{1 ≤ i ≤ s} b[i] * vi.
```
"""
mutable struct SequentialRungeKuttaNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  ui::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  ti::Real
  aii::Real
  J::Union{Nothing,AbstractMatrix}
  r::Union{Nothing,AbstractVector}
end

function Algebra.allocate_residual(
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, dt, ui = disop.ti, disop.dt, disop.ui
  !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
  usi = (ui, x)
  r = allocate_residual(disop.odeop, ti, usi, disop.odeopcache)
  !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
  r
end

function Algebra.residual!(
  r::AbstractVector,
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, dt, ui = disop.ti, disop.dt, disop.ui
  !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
  usi = (ui, x)
  residual!(r, disop.odeop, ti, usi, disop.odeopcache)
  !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
  r
end

function Algebra.allocate_jacobian(
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, dt, ui = disop.ti, disop.dt, disop.ui
  !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
  usi = (ui, x)
  J = allocate_jacobian(disop.odeop, ti, usi, disop.odeopcache)
  !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
  J
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, dt, ui = disop.ti, disop.dt, disop.ui
  ws = (disop.aii * dt, 1)
  !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
  usi = (ui, x)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, ti, usi, ws, disop.odeopcache)
  !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
  J
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, disop::SequentialRungeKuttaNonlinearOperator,
  disslvrcaches
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, dt, us0 = disop.t0, disop.dt, disop.us0
  ui, vs, tableau = disop.ui, disop.vs, disop.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  num_stages = length(c)
  is_masslinear = ODEOperatorType(odeop) <: AbstractMassLinearODE

  u0, = us0

  # Solve stages
  for i in 1:num_stages
    ti = t0 + c[i] * dt
    update_odeopcache!(odeopcache, odeop, ti)

    # Take linear combination of previous stages
    copy!(ui, u0)
    for j in 1:i-1
      coef = A[i, j]
      if !iszero(coef)
        axpy!(coef * dt, vs[j], ui)
      end
    end

    # Update operator state
    disop.ti = ti
    disop.aii = A[i, i]

    # Solve stage
    explicit = iszero(A[i, i])
    disslvr, islvr = get_solver_index(odeslvr, explicit)
    disslvrcache = disslvrcaches[islvr]

    vi = vs[i]
    fill!(vi, zero(eltype(vi)))
    if explicit && is_masslinear
      J, r = disop.J, disop.r

      usi = (ui, vi)
      w = 1
      filter = (true, true, false)

      fillstored!(J, zero(eltype(J)))
      jacobian!(J, odeop, ti, usi, 1, w, odeopcache)
      residual!(r, odeop, ti, usi, odeopcache; filter)
      rmul!(r, -1)

      _op = SequentialRungeKuttaLinearOperator(
        odeop, odeopcache,
        t0, us0, dt,
        ui, vs, tableau, J, r
      )
    else
      _op = disop
    end

    disslvrcache = solve!(vi, disslvr, _op, disslvrcache)
    disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)
  end

  # Express usF in terms of the solution of the discrete ODE operator
  usF = _finalize_rk!(usF, us0, vs, dt, b)

  (usF, disslvrcaches)
end

###################
# Linear operator #
###################
"""
    struct SequentialRungeKuttaLinearOperator <: LinearDiscreteODEOperator end

Linear operator corresponding to a sequential Runge-Kutta (explicit or
diagonally implicit) scheme:
```math
residual(ti, ui, vi) = mass(ti, ui) vi + res(ti, ui) = 0,

ti = t_n + c[i] * dt
ui = u_n + dt * ∑_{1 ≤ j < i} A[i, j] * vj + dt * A[i, i] x = 0
vi = x,

u_(n+1) = u_n + dt * ∑_{1 ≤ i ≤ s} b[i] * vi.
```
"""
struct SequentialRungeKuttaLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  ui::AbstractVector
  vs::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::SequentialRungeKuttaLinearOperator) = disop.J

Algebra.get_vector(disop::SequentialRungeKuttaLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, disop::SequentialRungeKuttaLinearOperator,
  disslvrcaches
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  t0, dt, us0 = disop.t0, disop.dt, disop.us0
  ui, vs, tableau = disop.ui, disop.vs, disop.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  num_stages = length(c)

  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  u0, = us0

  # Solve stages
  for i in 1:num_stages
    ti = t0 + c[i] * dt
    update_odeopcache!(odeopcache, odeop, ti)

    # Take linear combination of previous stages
    copy!(ui, u0)
    for j in 1:i-1
      coef = A[i, j]
      if !iszero(coef)
        axpy!(coef * dt, vs[j], ui)
      end
    end

    # Update jacobian and residual
    vi = vs[i]
    fill!(vi, zero(eltype(vi)))
    usi = (ui, vi) # vi does not matter
    ws = (A[i, i] * dt, 1)
    filter = (true, true, false)

    fillstored!(J, zero(eltype(J)))
    jacobians!(J, odeop, ti, usi, ws, odeopcache)
    residual!(r, odeop, ti, usi, odeopcache; filter)
    rmul!(r, -1)

    # Solve the discrete ODE operator
    disslvrcache = solve!(vi, disslvr, disop, disslvrcache)
  end

  # Express usF in terms of the solution of the discrete ODE operator
  usF = _finalize_rk!(usF, us0, vs, dt, b)

  disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)
  (usF, disslvrcaches)
end

#########
# Utils #
#########
function _finalize_rk!(
  usF::NTuple{1,AbstractVector}, us0::NTuple{1,AbstractVector},
  vs::AbstractVector, dt::Real, b::AbstractVector
)
  u0, = us0
  uF, = usF
  copy!(uF, u0)
  for i in eachindex(b)
    coef = b[i]
    if !iszero(coef)
      axpy!(coef * dt, vs[i], uF)
    end
  end
  (uF,)
end
