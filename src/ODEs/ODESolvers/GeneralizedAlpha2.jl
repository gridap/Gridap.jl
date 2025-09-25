"""
    struct GeneralizedAlpha2 <: ODESolver

Generalized-α second-order ODE solver.
```
residual(tx, ux, vx, ax) = 0,

tx = αf * t_n + (1 - αf) * t_(n+1)
ux = αf * u_n + (1 - αf) * u_(n+1)
vx = αf * v_n + (1 - αf) * v_(n+1)
ax = αm * a_n + (1 - αm) * a_(n+1),

u_(n+1) = u_n + dt * v_n + dt^2 / 2 * ((1 - 2 * β) * a_n + 2 * β * x)
v_(n+1) = v_n + dt * ((1 - γ) * a_n + γ * x)
a_(n+1) = x.
```
"""
struct GeneralizedAlpha2 <: ODESolver
  sysslvr::NonlinearSolver
  dt::Real
  αf::Real
  αm::Real
  γ::Real
  β::Real
end

# Constructors
function GeneralizedAlpha2(sysslvr::NonlinearSolver, dt::Real, ρ∞::Real)
  ρ∞01 = clamp(ρ∞, 0, 1)
  if ρ∞01 != ρ∞
    msg = """
    The parameter ρ∞ of the generalized-α scheme must lie between zero and one.
    Setting ρ∞ to $(ρ∞01).
    """
    @warn msg
    ρ∞ = ρ∞01
  end

  αf = ρ∞ / (1 + ρ∞)
  αm = (2 * ρ∞ - 1) / (1 + ρ∞)
  γ = 1 / 2 - αm + αf
  β = (1 - αm + αf)^2 / 4
  GeneralizedAlpha2(sysslvr, dt, αf, αm, γ, β)
end

"""
    Newmark(sysslvr::NonlinearSolver, dt::Real, γ::Real, β::Real)
"""
function Newmark(sysslvr::NonlinearSolver, dt::Real, γ::Real, β::Real)
  γ01 = clamp(γ, 0, 1)
  if γ01 != γ
    msg = """
    The parameter γ of the Newmark scheme must lie between zero and one.
    Setting γ to $(γ01).
    """
    @warn msg
    γ = γ01
  end

  β01 = clamp(β, 0, 1)
  if β01 != β
    msg = """
    The parameter β of the Newmark scheme must lie between zero and one.
    Setting β to $(β01).
    """
    @warn msg
    β = β01
  end

  αf, αm = 0.0, 0.0
  GeneralizedAlpha2(sysslvr, dt, αf, αm, γ, β)
end

# Default allocate_odecache without acceleration
function allocate_odecache(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  u0, v0 = us0[1], us0[2]
  allocate_odecache(odeslvr, odeop, t0, (u0, v0, v0))
end

##################
# Nonlinear case #
##################
function allocate_odecache(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  t0::Real, us0::NTuple{3,AbstractVector}
)
  u0, v0, a0 = us0[1], us0[2], us0[3]
  us0N = (u0, v0, a0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  uα, vα, aα = copy(u0), copy(v0), copy(a0)

  sysslvrcache = nothing
  odeslvrcache = (uα, vα, aα, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_start(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  t0::Real, us0::NTuple{2,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0 = us0[1], us0[2]
  odeslvrcache, odeopcache = odecache
  uα, vα, aα, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr

  # Allocate state
  s0, s1, s2 = copy(u0), copy(v0), copy(v0)

  # Define scheme
  x = s2
  tx = t0
  usx(x) = (u0, v0, x)
  ws = (0, 0, 1)

  # Update ODE operator cache
  update_odeopcache!(odeopcache, odeop, tx)

  # Create and solve stage operator
  stageop = NonlinearStageOperator(
    odeop, odeopcache,
    tx, usx, ws
  )

  sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  state0 = (s0, s1, s2)

  # Pack outputs
  odeslvrcache = (uα, vα, aα, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (state0, odecache)
end

function ode_start(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  t0::Real, us0::NTuple{3,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0, a0 = us0[1], us0[2], us0[3]

  # Allocate state
  s0, s1, s2 = copy(u0), copy(v0), copy(a0)

  # Update state
  state0 = (s0, s1, s2)

  # Pack outputs
  (state0, odecache)
end

function ode_march!(
  stateF::NTuple{3,AbstractVector},
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  t0::Real, state0::NTuple{3,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0, a0 = state0[1], state0[2], state0[3]
  odeslvrcache, odeopcache = odecache
  uα, vα, aα, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, αf, αm, γ, β = odeslvr.dt, odeslvr.αf, odeslvr.αm, odeslvr.γ, odeslvr.β

  # Define scheme
  tx = t0 + (1 - αf) * dt
  x = stateF[3]
  function usx(x)
    copy!(uα, u0)
    axpy!((1 - αf) * dt, v0, uα)
    axpy!((1 - αf) * (1 - 2 * β) * dt^2 / 2, a0, uα)
    axpy!((1 - αf) * β * dt^2, x, uα)

    copy!(vα, v0)
    axpy!((1 - αf) * (1 - γ) * dt, a0, vα)
    axpy!((1 - αf) * γ * dt, x, vα)

    copy!(aα, a0)
    rmul!(aα, αm)
    axpy!(1 - αm, x, aα)

    (uα, vα, aα)
  end
  ws = ((1 - αf) * β * dt^2, (1 - αf) * γ * dt, 1 - αm)

  # Update ODE operator cache
  update_odeopcache!(odeopcache, odeop, tx)

  # Create and solve stage operator
  stageop = NonlinearStageOperator(
    odeop, odeopcache,
    tx, usx, ws
  )

  sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  tF = t0 + dt
  stateF = _update_alpha2!(stateF, state0, dt, x, γ, β)

  # Pack outputs
  odeslvrcache = (uα, vα, aα, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

###############
# Linear case #
###############
function allocate_odecache(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{3,AbstractVector}
)
  u0, v0, a0 = us0[1], us0[2], us0[3]
  us0N = (u0, v0, a0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  uα, vα, aα = zero(u0), zero(v0), zero(a0)

  constant_stiffness = is_form_constant(odeop, 0)
  constant_damping = is_form_constant(odeop, 1)
  constant_mass = is_form_constant(odeop, 2)
  reuse = (constant_stiffness && constant_damping && constant_mass)

  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  sysslvrcache = nothing
  odeslvrcache = (reuse, uα, vα, aα, J, r, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_start(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{2,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0 = us0[1], us0[2]
  odeslvrcache, odeopcache = odecache
  reuse, uα, vα, aα, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr

  # Allocate state
  s0, s1, s2 = copy(u0), copy(v0), copy(v0)

  # Define scheme
  # Set x to zero to split jacobian and residual
  x = s2
  fill!(x, zero(eltype(x)))
  tx = t0
  usx = (u0, v0, x)
  ws = (0, 0, 1)

  # Update ODE operator cache
  update_odeopcache!(odeopcache, odeop, tx)

  # Create and solve stage operator
  stageop = LinearStageOperator(
    odeop, odeopcache,
    tx, usx, ws,
    J, r, false, sysslvrcache
  )

  sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  state0 = (s0, s1, s2)

  # Pack outputs
  odeslvrcache = (reuse, uα, vα, aα, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (state0, odecache)
end

function ode_start(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{3,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0, a0 = us0[1], us0[2], us0[3]

  # Allocate state
  s0, s1, s2 = copy(u0), copy(v0), copy(a0)

  # Update state
  state0 = (s0, s1, s2)

  # Pack outputs
  (state0, odecache)
end

function ode_march!(
  stateF::NTuple{3,AbstractVector},
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, state0::NTuple{3,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0, a0 = state0[1], state0[2], state0[3]
  odeslvrcache, odeopcache = odecache
  reuse, uα, vα, aα, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, αf, αm, γ, β = odeslvr.dt, odeslvr.αf, odeslvr.αm, odeslvr.γ, odeslvr.β

  # Define scheme
  x = stateF[3]
  tx = t0 + (1 - αf) * dt
  copy!(uα, u0)
  axpy!((1 - αf) * dt, v0, uα)
  axpy!((1 - αf) * (1 - 2 * β) * dt^2 / 2, a0, uα)
  copy!(vα, v0)
  axpy!((1 - αf) * (1 - γ) * dt, a0, vα)
  copy!(aα, a0)
  rmul!(aα, αm)
  usx = (uα, vα, aα)
  ws = ((1 - αf) * β * dt^2, (1 - αf) * γ * dt, 1 - αm)

  # Update ODE operator cache
  update_odeopcache!(odeopcache, odeop, tx)

  # Solve the discrete ODE operator
  stageop = LinearStageOperator(
    odeop, odeopcache,
    tx, usx, ws,
    J, r, reuse, sysslvrcache
  )

  sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  tF = t0 + dt
  stateF = _update_alpha2!(stateF, state0, dt, x, γ, β)

  # Pack outputs
  odeslvrcache = (reuse, uα, vα, aα, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _update_alpha2!(
  stateF::NTuple{3,AbstractVector}, state0::NTuple{3,AbstractVector},
  dt::Real, x::AbstractVector, γ::Real, β::Real
)
  # uF = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
  # vF = v0 + dt * ((1 - γ) * a0 + γ * x)
  # We always have x === aF
  u0, v0, a0 = state0[1], state0[2], state0[3]
  uF, vF, aF = stateF[1], stateF[2], stateF[3]

  copy!(uF, u0)
  axpy!(dt, v0, uF)
  axpy!((1 - 2 * β) * dt^2 / 2, a0, uF)
  axpy!(β * dt^2, x, uF)

  copy!(vF, v0)
  axpy!((1 - γ) * dt, a0, vF)
  axpy!(γ * dt, x, vF)

  (uF, vF, aF)
end
