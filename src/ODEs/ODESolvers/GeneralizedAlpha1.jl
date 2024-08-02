"""
    struct GeneralizedAlpha1 <: ODESolver

Generalized-α first-order ODE solver.
```
residual(tx, ux, vx) = 0,

tx = (1 - αf) * t_n + αf * t_(n+1)
ux = (1 - αf) * u_n + αf * u_(n+1)
vx = (1 - αm) * v_n + αm * v_(n+1),

u_(n+1) = u_n + dt * ((1 - γ) * v_n + γ * x)
v_(n+1) = x.
```
"""
struct GeneralizedAlpha1 <: ODESolver
  sysslvr::NonlinearSolver
  dt::Real
  αf::Real
  αm::Real
  γ::Real
end

# Constructors
function GeneralizedAlpha1(
  sysslvr::NonlinearSolver,
  dt::Real, ρ∞::Real
)
  ρ∞01 = clamp(ρ∞, 0, 1)
  if ρ∞01 != ρ∞
    msg = """
    The parameter ρ∞ of the generalized-α scheme must lie between zero and one.
    Setting ρ∞ to $(ρ∞01).
    """
    @warn msg
    ρ∞ = ρ∞01
  end

  αf = 1 / (1 + ρ∞)
  αm = (3 - ρ∞) / (1 + ρ∞) / 2
  γ = 1 / 2 + αm - αf

  GeneralizedAlpha1(sysslvr, dt, αf, αm, γ)
end

# Default allocate_odecache without velocity
function allocate_odecache(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  allocate_odecache(odeslvr, odeop, t0, (u0, u0))
end

##################
# Nonlinear case #
##################
function allocate_odecache(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  u0, v0 = us0[1], us0[2]
  us0N = (u0, v0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  uα, vα = copy(u0), copy(v0)

  sysslvrcache = nothing
  odeslvrcache = (uα, vα, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_start(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = us0[1]
  odeslvrcache, odeopcache = odecache
  uα, vα, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr

  # Allocate state
  s0, s1 = copy(u0), copy(u0)

  # Define scheme
  x = s1
  tx = t0
  usx(x) = (u0, x)
  ws = (0, 1)

  # Update ODE operator cache
  update_odeopcache!(odeopcache, odeop, tx)

  # Create and solve stage operator
  stageop = NonlinearStageOperator(
    odeop, odeopcache,
    tx, usx, ws
  )

  sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  state0 = (s0, s1)

  # Pack outputs
  odeslvrcache = (uα, vα, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (state0, odecache)
end

function ode_start(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  t0::Real, us0::NTuple{2,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0 = us0[1], us0[2]

  # Allocate state
  s0, s1 = copy(u0), copy(v0)

  # Update state
  state0 = (s0, s1)

  # Pack outputs
  (state0, odecache)
end

function ode_march!(
  stateF::NTuple{2,AbstractVector},
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  t0::Real, state0::NTuple{2,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0 = state0[1], state0[2]
  odeslvrcache, odeopcache = odecache
  uα, vα, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, αf, αm, γ = odeslvr.dt, odeslvr.αf, odeslvr.αm, odeslvr.γ

  # Define scheme
  x = stateF[2]
  tx = t0 + αf * dt
  function usx(x)
    # uα = u0 + αf * dt * [(1 - γ) * v0 + γ * x]
    copy!(uα, u0)
    axpy!(αf * (1 - γ) * dt, v0, uα)
    axpy!(αf * γ * dt, x, uα)

    # vα = (1 - αm) * v0 + αm * x
    copy!(vα, v0)
    rmul!(vα, 1 - αm)
    axpy!(αm, x, vα)

    (uα, vα)
  end
  ws = (αf * γ * dt, αm)

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
  stateF = _update_alpha1!(stateF, state0, dt, x, γ)

  # Pack outputs
  odeslvrcache = (uα, vα, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

###############
# Linear case #
###############
function allocate_odecache(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{2,AbstractVector}
)
  u0, v0 = us0[1], us0[2]
  us0N = (u0, v0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  uα, vα = zero(u0), zero(v0)

  constant_stiffness = is_form_constant(odeop, 0)
  constant_mass = is_form_constant(odeop, 1)
  reuse = (constant_stiffness && constant_mass)

  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  sysslvrcache = nothing
  odeslvrcache = (reuse, uα, vα, J, r, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_start(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = us0[1]
  odeslvrcache, odeopcache = odecache
  reuse, uα, vα, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr

  # Allocate state
  s0, s1 = copy(u0), copy(u0)

  # Define scheme
  # Set x to zero to split jacobian and residual
  x = s1
  fill!(x, zero(eltype(x)))
  tx = t0
  usx = (u0, x)
  ws = (0, 1)

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
  state0 = (s0, s1)

  # Pack outputs
  odeslvrcache = (reuse, uα, vα, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (state0, odecache)
end

function ode_start(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{2,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0 = us0[1], us0[2]

  # Allocate state
  s0, s1 = copy(u0), copy(v0)

  # Update state
  state0 = (s0, s1)

  # Pack outputs
  (state0, odecache)
end

function ode_march!(
  stateF::NTuple{2,AbstractVector},
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, state0::NTuple{2,AbstractVector},
  odecache
)
  # Unpack inputs
  u0, v0 = state0[1], state0[2]
  odeslvrcache, odeopcache = odecache
  reuse, uα, vα, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, αf, αm, γ = odeslvr.dt, odeslvr.αf, odeslvr.αm, odeslvr.γ

  # Define scheme
  dtα = αf * dt
  tx = t0 + dtα
  x = stateF[2]
  copy!(uα, u0)
  axpy!(αf * (1 - γ) * dt, v0, uα)
  copy!(vα, v0)
  rmul!(vα, 1 - αm)
  usx = (uα, vα)
  ws = (αf * γ * dt, αm)

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
  stateF = _update_alpha1!(stateF, state0, dt, x, γ)

  # Pack outputs
  odeslvrcache = (reuse, uα, vα, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _update_alpha1!(
  stateF::NTuple{2,AbstractVector}, state0::NTuple{2,AbstractVector},
  dt::Real, x::AbstractVector, γ::Real
)
  # uF = u0 + dt * ((1 - γ) * v0 + γ * x)
  # vF = x
  # We always have x === vF
  u0, v0 = state0[1], state0[2]
  uF, vF = stateF[1], stateF[2]

  copy!(uF, u0)
  axpy!((1 - γ) * dt, v0, uF)
  axpy!(γ * dt, x, uF)

  (uF, vF)
end
