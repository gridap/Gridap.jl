"""
    struct ThetaMethod <: ODESolver end

θ-method ODE solver.
"""
struct ThetaMethod <: ODESolver
  disslvr::NonlinearSolver
  dt::Real
  θ::Real

  function ThetaMethod(disslvr, dt, θ)
    θ01 = clamp(θ, 0, 1)
    if θ01 != θ
      msg = """
      The parameter θ of the θ-method must lie between zero and one.
      Setting θ to $(θ01).
      """
      @warn msg
    end

    if iszero(θ01)
      ForwardEuler(disslvr, dt)
    else
      new(disslvr, dt, θ01)
    end
  end
end

MidPoint(disslvr, dt) = ThetaMethod(disslvr, dt, 0.5)
BackwardEuler(disslvr, dt) = ThetaMethod(disslvr, dt, 1)

# ODESolver interface
function get_dt(odeslvr::ThetaMethod)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::ThetaMethod,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  (similar(x),)
end

function allocate_disopcache(
  odeslvr::ThetaMethod,
  odeop::ODEOperator{LinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  J = allocate_jacobian(odeop, t, (x, x), odeopcache)
  r = allocate_residual(odeop, t, (x, x), odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::ThetaMethod, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  tθ::Real, dtθ::Real
)
  usθ = disopcache
  ThetaMethodNonlinearOperator(
    odeop, odeopcache,
    tθ, us0, dt,
    dtθ, usθ
  )
end

function DiscreteODEOperator(
  odeslvr::ThetaMethod, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  tθ::Real, dtθ::Real
)
  J, r = disopcache
  ThetaMethodLinearOperator(
    odeop, odeopcache,
    tθ, us0, dt,
    dtθ, J, r
  )
end

###############
# solve_step! #
###############
function solve_step!(
  usF::NTuple{1,AbstractVector},
  odeslvr::ThetaMethod, odeop::ODEOperator,
  us0::NTuple{1,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, = us0
  dt = get_dt(odeslvr)
  θ = odeslvr.θ
  dtθ = θ * dt
  tθ = t0 + dtθ

  # Allocate or unpack cache
  if isnothing(cache)
    odeopcache = allocate_odeopcache(odeop, t0, (u0, u0))
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, u0)
    disslvrcache = allocate_disslvrcache(odeslvr)
  else
    odeopcache, disopcache, disslvrcache = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop,
    odeopcache, disopcache,
    t0, us0, dt,
    tθ, dtθ
  )

  # Solve discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (usF, tF, cache)
end

"""
    struct ThetaMethodNonlinearOperator <: DiscreteODEOperator

Nonlinear discrete operator corresponding to the θ-method scheme:
```math
residual(tθ, uθ, vθ) = 0,

tθ = t_n + θ * dt
uθ = u_n + θ * dt * x
vθ = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ThetaMethodNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  tθ::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  dtθ::Real
  usθ::NTuple{1,AbstractVector}
end

function Algebra.allocate_residual(
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  tθ, usθ, us0 = disop.tθ, disop.usθ, disop.us0
  dtθ = disop.dtθ
  usθ = _fill_usθ!(usθ, us0, x, dtθ)
  allocate_residual(disop.odeop, tθ, usθ, disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  tθ, usθ, us0 = disop.tθ, disop.usθ, disop.us0
  dtθ = disop.dtθ
  usθ = _fill_usθ!(usθ, us0, x, dtθ)
  residual!(r, disop.odeop, tθ, usθ, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  tθ, usθ, us0 = disop.tθ, disop.usθ, disop.us0
  dtθ = disop.dtθ
  usθ = _fill_usθ!(usθ, us0, x, dtθ)
  allocate_jacobian(disop.odeop, tθ, usθ, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  tθ, usθ, us0 = disop.tθ, disop.usθ, disop.us0
  dtθ = disop.dtθ
  usθ = _fill_usθ!(usθ, us0, x, dtθ)
  ws = _get_ws(dtθ)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tθ, usθ, ws, disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ThetaMethodNonlinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  tθ, dt, us0 = disop.tθ, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, tθ)

  vF, = usF
  fill!(vF, zero(eltype(vF)))

  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  _fill_usF!(usF, us0, vF, dt)

  (usF, disslvrcache)
end

"""
    struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the θ-method scheme:
```math
residual(tθ, uθ, vθ) = mass(tθ) vθ + stiffness(tθ) uθ + res(tθ) = 0,

tθ = t_n + θ * dt
uθ = u_n + θ * dt * x
vθ = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  tθ::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  dtθ::Real
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::ThetaMethodLinearOperator) = disop.J

Algebra.get_vector(disop::ThetaMethodLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ThetaMethodLinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  tθ, dt, us0 = disop.tθ, disop.dt, disop.us0
  dtθ = disop.dtθ

  update_odeopcache!(odeopcache, odeop, tθ)

  u0, = us0
  usθ = (u0, u0)
  ws = _get_ws(dtθ)

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, odeop, tθ, usθ, ws, odeopcache)
  residual!(r, odeop, tθ, usθ, odeopcache, include_highest=false)
  rmul!(r, -1)

  vF = usF[1]
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  usF = _fill_usF!(usF, us0, vF, dt)

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _fill_usF!(
  usF::NTuple{1,AbstractVector}, us0::NTuple{1,AbstractVector}, x,
  dt
)
  u0, = us0
  uF, = usF
  @. uF = u0 + dt * x
  (uF,)
end

function _get_ws(dtθ)
  wu = dtθ
  wv = 1
  (wu, wv)
end

function _fill_usθ!(
  usθ::NTuple{1,AbstractVector}, us0::NTuple{1,AbstractVector}, x,
  dtθ
)
  u0, = us0
  uθ, = usθ
  @. uθ = u0 + dtθ * x
  vθ = x
  (uθ, vθ)
end
