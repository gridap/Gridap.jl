"""
    struct ForwardEuler <: ODESolver end

Forward Euler ODE solver.
```
residual(tx, ux, vx) = 0,

tx = t_n
ux = u_n
vx = x,

u_(n+1) = u_n + dt * x.
```
"""
struct ForwardEuler <: ODESolver
  sysslvr::NonlinearSolver
  dt::Real
end

##################
# Nonlinear case #
##################
function allocate_odecache(
  odeslvr::ForwardEuler, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  sysslvrcache = nothing
  odeslvrcache = (sysslvrcache,)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::ForwardEuler, odeop::ODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  sysslvrcache, = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt = odeslvr.dt

  # Define scheme
  x = stateF[1]
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
  tF = t0 + dt
  stateF = _update_euler!(stateF, state0, dt, x)

  # Pack outputs
  odeslvrcache = (sysslvrcache,)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

###############
# Linear case #
###############
function allocate_odecache(
  odeslvr::ForwardEuler, odeop::ODEOperator{<:AbstractQuasilinearODE},
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  is_semilinear = (ODEOperatorType(odeop) <: AbstractSemilinearODE)
  constant_mass = is_form_constant(odeop, 1)
  reuse = (is_semilinear && constant_mass)

  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  sysslvrcache = nothing
  odeslvrcache = (reuse, J, r, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::ForwardEuler, odeop::ODEOperator{<:AbstractQuasilinearODE},
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  reuse, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt = odeslvr.dt

  # Define scheme
  # Set x to zero to split jacobian and residual
  x = stateF[1]
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
    J, r, reuse, sysslvrcache
  )

  sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  tF = t0 + dt
  stateF = _update_euler!(stateF, state0, dt, x)

  # Pack outputs
  odeslvrcache = (reuse, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _update_euler!(
  stateF::NTuple{1,AbstractVector}, state0::NTuple{1,AbstractVector},
  dt::Real, x::AbstractVector
)
  # uF = u0 + dt * x
  # We always have x === uF
  u0 = state0[1]
  uF = stateF[1]

  rmul!(uF, dt)
  axpy!(1, u0, uF)

  (uF,)
end
