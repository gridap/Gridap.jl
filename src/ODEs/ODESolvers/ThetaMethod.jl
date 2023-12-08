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
  t0::Real, us0::NTuple{2,AbstractVector}
)
  zero(us0[1])
end

function allocate_disopcache(
  odeslvr::ThetaMethod,
  odeop::ODEOperator{<:AbstractLinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::ThetaMethod, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real, tθ::Real, dtθ::Real
)
  uθ = disopcache
  ThetaMethodNonlinearOperator(
    odeop, odeopcache,
    u0, dt, tθ, uθ, dtθ
  )
end

function DiscreteODEOperator(
  odeslvr::ThetaMethod, odeop::ODEOperator{<:AbstractLinearODE},
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real, tθ::Real, dtθ::Real
)
  J, r = disopcache
  ThetaMethodLinearOperator(
    odeop, odeopcache,
    u0, dt, tθ, dtθ,
    J, r
  )
end

function solve_odeop!(
  usF::NTuple{1,AbstractVector},
  odeslvr::ThetaMethod, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector},
  cache
)
  u0 = us0[1]
  dt = get_dt(odeslvr)
  θ = odeslvr.θ
  dtθ = θ * dt
  tθ = t0 + dtθ
  tF = t0 + dt

  # Allocate or unpack cache
  if isnothing(cache)
    us0_full = (u0, u0)
    odeopcache = allocate_odeopcache(odeop, t0, us0_full)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, us0_full)
    disslvrcache = allocate_disslvrcache(odeslvr)
  else
    odeopcache, disopcache, disslvrcache = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop, odeopcache, disopcache,
    t0, u0, dt, tθ, dtθ
  )

  # Solve the discrete ODE operator
  usF, disslvrcache = solve_disop!(usF, odeslvr.disslvr, disop, disslvrcache)
  cache = (odeopcache, disopcache, disslvrcache)

  (tF, usF, cache)
end

######################
# Nonlinear operator #
######################
"""
    struct ThetaMethodNonlinearOperator <: DiscreteODEOperator

Nonlinear discrete operator corresponding to the θ-method scheme:
```math
residual(tx, ux, vx) = 0,

tx = t_n + θ * dt
ux = u_n + θ * dt * x
vx = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ThetaMethodNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  u0::AbstractVector
  dt::Real
  tθ::Real
  uθ::AbstractVector
  dtθ::Real
end

function Algebra.residual!(
  r::AbstractVector,
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  u0 = disop.u0
  tθ, uθ, dtθ = disop.tθ, disop.uθ, disop.dtθ
  # Residual: (u0 + dtθ * x, x)
  copy!(uθ, u0)
  axpy!(dtθ, x, uθ)

  tx = tθ
  usx = (uθ, x)
  residual!(r, disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  tθ, uθ = disop.tθ, disop.uθ
  tx = tθ
  usx = (uθ, x)
  allocate_jacobian(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::ThetaMethodNonlinearOperator,
  x::AbstractVector
)
  u0 = disop.u0
  tθ, uθ, dtθ = disop.tθ, disop.uθ, disop.dtθ
  # Jacobian: (u0 + dtθ * x, x)
  copy!(uθ, u0)
  axpy!(dtθ, x, uθ)

  tx = tθ
  usx = (uθ, x)
  ws = (dtθ, 1)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tx, usx, ws, disop.odeopcache)
end

function solve_disop!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ThetaMethodNonlinearOperator,
  disslvrcache
)
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  u0, dt = disop.u0, disop.dt
  tθ = disop.tθ

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, tθ)

  # Solve the discrete ODE operator
  x = uF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  uF = _finalize_theta!(uF, u0, dt)
  usF = (uF,)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the θ-method scheme:
```math
residual(tx, ux, vx) = mass(tx) vx + stiffness(tx) ux + res(tx) = 0,

tx = t_n + θ * dt
ux = u_n + θ * dt * x
vx = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ThetaMethodLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  u0::AbstractVector
  dt::Real
  tθ::Real
  dtθ::Real
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::ThetaMethodLinearOperator) = disop.J

Algebra.get_vector(disop::ThetaMethodLinearOperator) = disop.r

# The jacobian matrix is constant if the ODE operator is linear and has
# constant forms
function is_jacobian_constant(disop::ThetaMethodLinearOperator)
  odeop = disop.odeop
  constant_jacobian = false
  if ODEOperatorType(odeop) <: AbstractLinearODE
    constant_stiffness = is_form_constant(odeop, 0)
    constant_mass = is_form_constant(odeop, 1)
    constant_jacobian = constant_stiffness && constant_mass
  end
  constant_jacobian
end

function solve_disop!(
  usF::NTuple{1,AbstractVector},
  disslvr::NonlinearSolver, disop::ThetaMethodLinearOperator,
  disslvrcache
)
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  u0, dt = disop.u0, disop.dt
  tθ, dtθ = disop.tθ, disop.dtθ

  J, r = disop.J, disop.r

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, tθ)

  # Residual: (u0 + dtθ * x, x)
  # Jacobian: (u0 + dtθ * x, x)
  # Take x = 0 to split the mass term from the residual
  # Trick: use uF to store x = 0
  x = uF
  fill!(x, zero(eltype(x)))

  tx = tθ
  usx = (u0, x)
  ws = (dtθ, 1)
  # If the jacobian is constant, the following call will only retrieve the
  # jacobian from the ODEOpFromFEOpCache.
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, odeop, tx, usx, ws, odeopcache)
  residual!(r, odeop, tx, usx, odeopcache)
  rmul!(r, -1)

  # Solve the discrete ODE operator
  x = uF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalize
  uF = _finalize_theta!(uF, u0, dt)
  usF = (uF,)

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _finalize_theta!(uF::AbstractVector, u0::AbstractVector, dt::Real)
  # @. uF = u0 + dt * x
  axpby!(1, u0, dt, uF)
  uF
end
