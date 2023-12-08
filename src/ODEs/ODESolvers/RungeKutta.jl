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

function solve_odeop!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector},
  cache
)
  u0 = us0[1]
  dt = get_dt(odeslvr)
  tF = t0 + dt
  tableau = get_tableau(odeslvr)
  num_stages = length(get_nodes(tableau))

  # Allocate or unpack cache
  if isnothing(cache)
    us0_full = (u0, u0)
    odeopcache = allocate_odeopcache(odeop, t0, us0_full)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, us0_full)
    disslvrcache = allocate_disslvrcache(odeslvr)
    slopes = [zero(u0) for _ in 1:num_stages]
  else
    odeopcache, disopcache, disslvrcache, slopes = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop, odeopcache, disopcache,
    t0, u0, dt, slopes, tableau
  )

  # Solve the discrete ODE operator
  usF, disslvrcache = solve_disop!(usF, odeslvr, disop, disslvrcache)
  cache = (odeopcache, disopcache, disslvrcache, slopes)

  (tF, usF, cache)
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
  t0::Real, us0::NTuple{2,AbstractVector}
)
  ui = zero(us0[1])
  J, r = nothing, nothing
  (ui, J, r)
end

function allocate_disopcache(
  odeslvr::EXRungeKutta,
  odeop::ODEOperator{<:AbstractQuasilinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  ui = zero(us0[1])
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (ui, J, r)
end

function allocate_disslvrcache(odeslvr::EXRungeKutta)
  (nothing,)
end

function DiscreteODEOperator(
  odeslvr::EXRungeKutta, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
  slopes::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  ti, aii = t0, zero(t0)
  SequentialRungeKuttaNonlinearOperator(
    odeop, odeopcache,
    t0, u0, dt, ti, ui, aii,
    slopes, tableau, J, r
  )
end

function DiscreteODEOperator(
  odeslvr::EXRungeKutta, odeop::ODEOperator{<:AbstractQuasilinearODE},
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
  slopes::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  SequentialRungeKuttaLinearOperator(
    odeop, odeopcache,
    t0, u0, dt, ui,
    slopes, tableau, J, r,
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
  t0::Real, us0::NTuple{2,AbstractVector}
)
  ui = zero(us0[1])
  J, r = nothing, nothing
  (ui, J, r)
end

function allocate_disopcache(
  odeslvr::DIMRungeKutta,
  odeop::ODEOperator{<:AbstractQuasilinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  ui = zero(us0[1])
  tableau = get_tableau(odeslvr)
  A = get_matrix(tableau)
  if any(i -> iszero(A[i, i]), axes(A, 2))
    J = allocate_jacobian(odeop, t0, us0, odeopcache)
    r = allocate_residual(odeop, t0, us0, odeopcache)
  else
    J, r = nothing, nothing
  end
  (ui, J, r)
end

function allocate_disopcache(
  odeslvr::DIMRungeKutta,
  odeop::ODEOperator{<:AbstractLinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  ui = zero(us0[1])
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (ui, J, r)
end

function allocate_disslvrcache(odeslvr::DIMRungeKutta)
  (nothing, nothing)
end

function DiscreteODEOperator(
  odeslvr::DIMRungeKutta, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
  slopes::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  ti, aii = t0, zero(t0)
  SequentialRungeKuttaNonlinearOperator(
    odeop, odeopcache,
    t0, u0, dt, ti, ui, aii,
    slopes, tableau, J, r
  )
end

function DiscreteODEOperator(
  odeslvr::DIMRungeKutta, odeop::ODEOperator{<:AbstractLinearODE},
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
  slopes::AbstractVector{<:AbstractVector}, tableau::AbstractTableau
)
  ui, J, r = disopcache
  SequentialRungeKuttaLinearOperator(
    odeop, odeopcache,
    t0, u0, dt, ui,
    slopes, tableau, J, r
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
residual(tx, ux, vx) = 0,

tx = t_n + c[i] * dt
ux = u_n + dt * ∑_{1 ≤ j < i} A[i, j] * slopes[j] + dt * A[i, i] * dt * x
vx = x,

u_(n+1) = u_n + dt * ∑_{1 ≤ i ≤ s} b[i] * slopes[i].
```
"""
mutable struct SequentialRungeKuttaNonlinearOperator{C} <: DiscreteODEOperator
  odeop::ODEOperator{C}
  odeopcache
  t0::Real
  u0::AbstractVector
  dt::Real
  ti::Real
  ui::AbstractVector
  aii::Real
  slopes::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::Union{Nothing,AbstractMatrix}
  r::Union{Nothing,AbstractVector}
end

function Algebra.residual!(
  r::AbstractVector,
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  dt = disop.dt
  ti, ui, aii = disop.ti, disop.ui, disop.aii
  # Residual: (ui + aii * dt * x, x)
  axpy!(aii * dt, x, ui)

  tx = ti
  usx = (ui, x)
  residual!(r, disop.odeop, tx, usx, disop.odeopcache)
  axpy!(-aii * dt, x, ui)

  r
end

function Algebra.allocate_jacobian(
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, ui = disop.ti, disop.ui
  tx = ti
  usx = (ui, x)
  allocate_jacobian(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::SequentialRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  dt = disop.dt
  ti, ui, aii = disop.ti, disop.ui, disop.aii
  # Residual: (ui + aii * dt * x, x)
  axpy!(aii * dt, x, ui)

  tx = ti
  usx = (ui, x)
  ws = (aii * dt, 1)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tx, usx, ws, disop.odeopcache)
  axpy!(-aii * dt, x, ui)
  J
end

function solve_disop!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, disop::SequentialRungeKuttaNonlinearOperator{C},
  disslvrcaches
) where {C}
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, u0, dt = disop.t0, disop.u0, disop.dt
  ui = disop.ui

  slopes, tableau = disop.slopes, disop.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  is_quasilinear = C <: AbstractQuasilinearODE

  # Solve stages
  for i in eachindex(c)
    ti = t0 + c[i] * dt
    aii = A[i, i]

    # Update the cache of the ODE operator (typically Dirichlet BCs)
    update_odeopcache!(odeopcache, odeop, ti)

    # Update solver state (ti, aii and ui)
    disop.ti = ti
    disop.aii = aii
    copy!(ui, u0)
    for j in 1:i-1
      coef = A[i, j]
      if !iszero(coef)
        axpy!(coef * dt, slopes[j], ui)
      end
    end

    # Decide whether the operator is linear or nonlinear
    explicit = iszero(aii)
    disslvr, islvr = get_solver_index(odeslvr, explicit)
    disslvrcache = disslvrcaches[islvr]

    if explicit && is_quasilinear
      J, r = disop.J, disop.r

      # Residual: (ui, x)
      # Jacobian: (ui, x)
      # Take x = 0 to split the mass term from the residual
      # Use slopes[i] to store x = 0
      x = slopes[i]
      fill!(x, zero(eltype(x)))

      tx = ti
      usx = (ui, x)
      w = 1
      fillstored!(J, zero(eltype(J)))
      jacobian!(J, odeop, tx, usx, 1, w, odeopcache)
      residual!(r, odeop, tx, usx, odeopcache)
      rmul!(r, -1)

      _op = SequentialRungeKuttaLinearOperator(
        odeop, odeopcache,
        t0, u0, dt, ui,
        slopes, tableau, J, r
      )
    else
      _op = disop
    end

    # Solve the discrete ODE operator
    x = slopes[i]
    disslvrcache = solve!(x, disslvr, _op, disslvrcache)
    disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)
  end

  # Finalize
  uF = _finalize_rk!(uF, u0, dt, slopes, b)
  usF = (uF,)

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
residual(tx, ux, vx) = mass(tx, ux) vx + res(tx, ux) = 0,

tx = t_n + c[i] * dt
ux = u_n + dt * ∑_{1 ≤ j < i} A[i, j] * slopes[j] + dt * A[i, i] x = 0
vx = x,

u_(n+1) = u_n + dt * ∑_{1 ≤ i ≤ s} b[i] * slopes[i].
```
"""
struct SequentialRungeKuttaLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  u0::AbstractVector
  dt::Real
  ui::AbstractVector
  slopes::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::SequentialRungeKuttaLinearOperator) = disop.J

Algebra.get_vector(disop::SequentialRungeKuttaLinearOperator) = disop.r

# For now only indicate that the jacobian matrix is constant when the
# tableau is explicit and the ODE operator is semilinear with constant mass.
# TODO as indicated in ODESolvers.jl, we should also save the jacobian matrices
# when the ODE operator is linear with constant forms.
function is_jacobian_constant(disop::SequentialRungeKuttaLinearOperator)
  odeop = disop.odeop
  constant_jacobian = false
  explicit_tableau = TableauType(disop.tableau) <: ExplicitTableau
  semilinear_ode = ODEOperatorType(odeop) <: AbstractSemilinearODE
  if explicit_tableau && semilinear_ode
    constant_jacobian = is_form_constant(odeop, 1)
  end
  constant_jacobian
end

function solve_disop!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, disop::SequentialRungeKuttaLinearOperator,
  disslvrcaches
)
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, u0, dt = disop.t0, disop.u0, disop.dt
  ui = disop.ui
  slopes, tableau = disop.slopes, disop.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)
  J, r = disop.J, disop.r

  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  # Solve stages
  for i in eachindex(c)
    ti = t0 + c[i] * dt

    # Update the cache of the ODE operator (typically Dirichlet BCs)
    update_odeopcache!(odeopcache, odeop, ti)

    # Residual: (ui, x)
    # Jacobian: (ui, x)
    # Take x = 0 to split the mass term from the residual
    # Use slopes[i] to store x = 0
    copy!(ui, u0)
    for j in 1:i-1
      coef = A[i, j]
      if !iszero(coef)
        axpy!(coef * dt, slopes[j], ui)
      end
    end
    x = slopes[i]
    fill!(x, zero(eltype(x)))

    tx = ti
    usx = (ui, x)
    ws = (A[i, i] * dt, 1)
    # If the jacobian is constant, the following call will only retrieve the
    # jacobian from the ODEOpFromFEOpCache.
    fillstored!(J, zero(eltype(J)))
    jacobians!(J, odeop, tx, usx, ws, odeopcache)
    residual!(r, odeop, tx, usx, odeopcache)
    rmul!(r, -1)

    # Solve the discrete ODE operator
    x = slopes[i]
    disslvrcache = solve!(x, disslvr, disop, disslvrcache)
  end
  disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)

  # Finalize
  uF = _finalize_rk!(uF, u0, dt, slopes, b)
  usF = (uF,)

  (usF, disslvrcaches)
end

#############
# Finalizer #
#############
function _finalize_rk!(
  uF::AbstractVector, u0::AbstractVector, dt::Real,
  slopes::AbstractVector, b::AbstractVector
)
  # @. uF = u0 + ∑_{1 ≤ i ≤ s} b[i] * dt * slopes[i]
  copy!(uF, u0)
  for i in eachindex(b)
    coef = b[i]
    if !iszero(coef)
      axpy!(coef * dt, slopes[i], uF)
    end
  end
  uF
end
