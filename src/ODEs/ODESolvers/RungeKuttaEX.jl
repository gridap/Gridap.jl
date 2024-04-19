################
# EXRungeKutta #
################
"""
    struct EXRungeKutta <: ODESolver end

Explicit Runge-Kutta ODE solver.
```
residual(tx, ux, vx) = 0,

tx = t_n + c[i] * dt
ux = u_n + ∑_{1 ≤ j < i} A[i, j] * dt * slopes[j]
vx = x
slopes[i] = x,

u_(n+1) = u_n + ∑_{1 ≤ i ≤ s} b[i] * dt * slopes[i].
```
"""
struct EXRungeKutta <: ODESolver
  sysslvr::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{ExplicitTableau}
end

##################
# Nonlinear case #
##################
function allocate_odecache(
  odeslvr::EXRungeKutta, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  ui_pre = zero(u0)
  num_stages = length(get_nodes(odeslvr.tableau))
  slopes = [zero(u0) for _ in 1:num_stages]

  sysslvrcache = nothing
  odeslvrcache = (ui_pre, slopes, sysslvrcache)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::EXRungeKutta, odeop::ODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  ui_pre, slopes, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, tableau = odeslvr.dt, odeslvr.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)

  for i in eachindex(c)
    # Define scheme
    x = slopes[i]
    tx = t0 + c[i] * dt
    copy!(ui_pre, u0)
    for j in 1:i-1
      axpy!(A[i, j] * dt, slopes[j], ui_pre)
    end
    usx(x) = (ui_pre, x)
    ws = (0, 1)

    # Update ODE operator cache
    update_odeopcache!(odeopcache, odeop, tx)

    # Create and solve stage operator
    stageop = NonlinearStageOperator(
      odeop, odeopcache,
      tx, usx, ws
    )

    sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)
  end

  # Update state
  tF = t0 + dt
  stateF = _update_exrk!(stateF, state0, dt, slopes, b)

  # Pack outputs
  odeslvrcache = (ui_pre, slopes, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

###############
# Linear case #
###############
function allocate_odecache(
  odeslvr::EXRungeKutta, odeop::ODEOperator{<:AbstractQuasilinearODE},
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  ui_pre = zero(u0)
  num_stages = length(get_nodes(odeslvr.tableau))
  slopes = [zero(u0) for _ in 1:num_stages]

  is_semilinear = ODEOperatorType(odeop) <: AbstractSemilinearODE
  constant_mass = is_form_constant(odeop, 1)
  reuse = (is_semilinear && constant_mass)

  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  sysslvrcache = nothing
  odeslvrcache = (reuse, ui_pre, slopes, J, r, sysslvrcache)

  (odeslvrcache, odeopcache)
end


function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::EXRungeKutta, odeop::ODEOperator{<:AbstractQuasilinearODE},
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  reuse, ui_pre, slopes, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, tableau = odeslvr.dt, odeslvr.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)

  for i in eachindex(c)
    # Define scheme
    # Set x to zero to split jacobian and residual
    x = slopes[i]
    fill!(x, zero(eltype(x)))
    tx = t0 + c[i] * dt
    copy!(ui_pre, u0)
    for j in 1:i-1
      axpy!(A[i, j] * dt, slopes[j], ui_pre)
    end
    usx = (ui_pre, x)
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
  end

  # Update state
  tF = t0 + dt
  stateF = _update_exrk!(stateF, state0, dt, slopes, b)

  # Pack outputs
  odeslvrcache = (reuse, ui_pre, slopes, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _update_exrk!(
  stateF::NTuple{1,AbstractVector}, state0::NTuple{1,AbstractVector},
  dt::Real, slopes::AbstractVector, b::AbstractVector
)
  # uF = u0 + ∑_{1 ≤ i ≤ s} b[i] * dt * slopes[i]
  u0 = state0[1]
  uF = stateF[1]

  copy!(uF, u0)
  for (bi, slopei) in zip(b, slopes)
    axpy!(bi * dt, slopei, uF)
  end

  (uF,)
end
