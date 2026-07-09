"""
    struct ThetaMethod <: ODESolver end

θ-method ODE solver.
```
residual(tx, ux, vx) = 0,

tx = t_n + θ * dt
ux = u_n + θ * dt * x
vx = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ThetaMethod <: ODESolver
  sysslvr::NonlinearSolver
  dt::Real
  θ::Real

  function ThetaMethod(sysslvr, dt, θ)
    θ01 = clamp(θ, 0, 1)
    if θ01 != θ
      msg = """
      The parameter θ of the θ-method must lie between zero and one.
      Setting θ to $(θ01).
      """
      @warn msg
    end

    if iszero(θ01)
      ForwardEuler(sysslvr, dt)
    else
      new(sysslvr, dt, θ01)
    end
  end
end

"""
    MidPoint(sysslvr, dt) = ThetaMethod(sysslvr, dt, 0.5)
"""
MidPoint(sysslvr, dt) = ThetaMethod(sysslvr, dt, 0.5)
"""
    BackwardEuler(sysslvr, dt) = ThetaMethod(sysslvr, dt, 1)
"""
BackwardEuler(sysslvr, dt) = ThetaMethod(sysslvr, dt, 1)

##################
# Nonlinear case #
##################
function allocate_odecache(
  odeslvr::ThetaMethod, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  uθ = copy(u0)

  sysslvrcache = nothing
  odeslvrcache = (uθ, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::ThetaMethod, odeop::ODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  uθ, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, θ = odeslvr.dt, odeslvr.θ

  # Define scheme
  x = stateF[1]
  dtθ = θ * dt
  tx = t0 + dtθ
  function usx(x)
    copy!(uθ, u0)
    axpy!(dtθ, x, uθ)
    (uθ, x)
  end
  ws = (dtθ, 1)

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
  stateF = _udate_theta!(stateF, state0, dt, x)

  # Pack outputs
  odeslvrcache = (uθ, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

###############
# Linear case #
###############
function allocate_odecache(
  odeslvr::ThetaMethod, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  constant_stiffness = is_form_constant(odeop, 0)
  constant_mass = is_form_constant(odeop, 1)
  reuse = (constant_stiffness && constant_mass)

  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  sysslvrcache = nothing
  odeslvrcache = (reuse, J, r, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::ThetaMethod, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  reuse, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, θ = odeslvr.dt, odeslvr.θ

  # Define scheme
  # Set x to zero to split jacobian and residual
  x = stateF[1]
  fill!(x, zero(eltype(x)))
  dtθ = θ * dt
  tx = t0 + dtθ
  usx = (u0, x)
  ws = (dtθ, 1)

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
  stateF = _udate_theta!(stateF, state0, dt, x)

  # Pack outputs
  odeslvrcache = (reuse, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _udate_theta!(
  stateF::NTuple{1,AbstractVector}, state0::NTuple{1,AbstractVector},
  dt::Real, x::AbstractVector
)
  # uF = u0 + dt * x
  # We always have x === uF
  u0 = state0[1]
  uF = stateF[1]
  rmul!(uF, dt)
  axpy!(1, u0, uF)
  stateF = (uF,)
end
