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
  similar(x)
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
  ulc = disopcache
  ThetaMethodNonlinearOperator(odeop, odeopcache, ulc, t0, us0, dt, tθ, dtθ)
end

function DiscreteODEOperator(
  odeslvr::ThetaMethod, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  tθ::Real, dtθ::Real
)
  J, r = disopcache
  ThetaMethodLinearOperator(odeop, odeopcache, J, r, t0, us0, dt, tθ, dtθ)
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
residual(t_n + θ * dt, u_n + θ * dt * x, x) = 0,
u_n+1 = u_n + dt * x
```
"""
struct ThetaMethodNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  ulc::AbstractVector
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  tθ::Real
  dtθ::Real
end

function Algebra.allocate_residual(
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  t, u, (u0,) = disop.tθ, disop.ulc, disop.us0
  u = _u_from_v!(u, u0, disop.dtθ, x)
  allocate_residual(disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  t, u, (u0,) = disop.tθ, disop.ulc, disop.us0
  u = _u_from_v!(u, u0, disop.dtθ, x)
  residual!(r, disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  t, u, (u0,) = disop.tθ, disop.ulc, disop.us0
  u = _u_from_v!(u, u0, disop.dtθ, x)
  allocate_jacobian(disop.odeop, t, (u, x), disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  t, u, (u0,) = disop.tθ, disop.ulc, disop.us0
  u = _u_from_v!(u, u0, disop.dtθ, x)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, t, (u, x), (disop.dtθ, 1), disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ThetaMethodNonlinearOperator,
  disslvrcache
)
  uF, = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  tθ, dt, (u0,) = disop.tθ, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, tθ)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  _u_from_v!(uF, u0, dt, vF)

  usF = (uF,)
  (usF, disslvrcache)
end

"""
    struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the θ-method scheme:
```math
residual(t_n + θ * dt, u_n + θ * dt * x, x) = mass(t_n + θ * dt) x + stiffness(t_n + θ * dt) (u_n + θ * dt * x) + res(t_n + θ * dt) = 0,
u_n+1 = u_n + dt * x
```
"""
struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  tθ::Real
  dtθ::Real
end

Algebra.get_matrix(disop::ThetaMethodLinearOperator) = disop.J

Algebra.get_vector(disop::ThetaMethodLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ThetaMethodLinearOperator,
  disslvrcache
)
  uF, = usF
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  tθ, dt, (u0,), dtθ = disop.tθ, disop.dt, disop.us0, disop.dtθ

  update_odeopcache!(odeopcache, odeop, tθ)

  u, x = u0, u0
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, odeop, tθ, (u, x), (dtθ, 1), odeopcache)
  residual!(r, odeop, tθ, (u, x), odeopcache, include_highest=false)
  rmul!(r, -1)

  vF = uF
  fill!(vF, zero(eltype(vF)))
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  _u_from_v!(uF, u0, dt, vF)

  usF = (uF,)
  (usF, disslvrcache)
end
