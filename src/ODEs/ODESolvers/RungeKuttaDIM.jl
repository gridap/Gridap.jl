
#################
# DIMRungeKutta #
#################
"""
    struct DIMRungeKutta <: ODESolver end

Diagonally-implicit Runge-Kutta ODE solver.
```
residual(tx, ux, vx) = 0,

tx = t_n + c[i] * dt
ux = u_n + dt * ∑_{1 ≤ j < i} A[i, j] * slopes[j] + dt * A[i, i] * x
vx = x,

u_(n+1) = u_n + dt * ∑_{1 ≤ i ≤ s} b[i] * slopes[i].
```
"""
struct DIMRungeKutta <: ODESolver
  sysslvr_nl::NonlinearSolver
  sysslvr_l::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{DiagonallyImplicitTableau}
end

##################
# Nonlinear case #
##################
function allocate_odecache(
  odeslvr::DIMRungeKutta, odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  ui_pre, ui = zero(u0), zero(u0)
  num_stages = length(get_nodes(odeslvr.tableau))
  slopes = [zero(u0) for _ in 1:num_stages]

  odeoptype = ODEOperatorType(odeop)
  has_explicit = odeoptype <: AbstractQuasilinearODE
  is_semilinear = odeoptype <: AbstractSemilinearODE
  mass_constant = is_form_constant(odeop, 1)
  reuse = (is_semilinear && mass_constant)

  J, r = nothing, nothing
  if has_explicit
    # Allocate J, r if there are explicit stages
    A = get_matrix(odeslvr.tableau)
    if any(i -> iszero(A[i, i]), axes(A, 2))
      J = allocate_jacobian(odeop, t0, us0N, odeopcache)
      r = allocate_residual(odeop, t0, us0N, odeopcache)
    end
  end

  sysslvrcaches = (nothing, nothing)
  odeslvrcache = (reuse, has_explicit, ui_pre, ui, slopes, J, r, sysslvrcaches)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::DIMRungeKutta, odeop::ODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  reuse, has_explicit, ui_pre, ui, slopes, J, r, sysslvrcaches = odeslvrcache
  sysslvrcache_nl, sysslvrcache_l = sysslvrcaches

  # Unpack solver
  sysslvr_nl, sysslvr_l = odeslvr.sysslvr_nl, odeslvr.sysslvr_l
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

    # Update ODE operator cache
    update_odeopcache!(odeopcache, odeop, tx)

    # Decide whether the stage is explicit or implicit
    # The stage becomes explicit when aii = 0 and the operator is quasilinear,
    # which is precomputed in has_explicit
    aii = A[i, i]
    explicit_stage = iszero(aii) && has_explicit

    if explicit_stage
      # Define scheme
      # Set x to zero to split jacobian and residual
      fill!(x, zero(eltype(x)))
      usx = (ui_pre, x)
      ws = (0, 1)

      # Create and solve stage operator
      stageop = LinearStageOperator(
        odeop, odeopcache,
        tx, usx, ws,
        J, r, reuse, sysslvrcache_l
      )

      sysslvrcache_l = solve!(x, sysslvr_l, stageop, sysslvrcache_l)
    else
      # Define scheme
      function usx(x)
        copy!(ui, ui_pre)
        axpy!(aii * dt, x, ui)
        (ui, x)
      end
      ws = (aii * dt, 1)

      # Create and solve stage operator
      stageop = NonlinearStageOperator(
        odeop, odeopcache,
        tx, usx, ws
      )

      sysslvrcache_nl = solve!(x, sysslvr_nl, stageop, sysslvrcache_nl)
    end
  end

  # Update state
  tF = t0 + dt
  stateF = _update_dimrk!(stateF, state0, dt, slopes, b)

  # Pack outputs
  sysslvrcaches = (sysslvrcache_nl, sysslvrcache_l)
  odeslvrcache = (reuse, has_explicit, ui_pre, ui, slopes, J, r, sysslvrcaches)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

###############
# Linear case #
###############
function allocate_odecache(
  odeslvr::DIMRungeKutta, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  ui_pre = zero(u0)
  num_stages = length(get_nodes(odeslvr.tableau))
  slopes = [zero(u0) for _ in 1:num_stages]

  stiffness_constant = is_form_constant(odeop, 0)
  mass_constant = is_form_constant(odeop, 1)
  reuse = (stiffness_constant && mass_constant)

  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  # Numerical setups for the linear solver
  # * If the mass and stiffness matrices are constant, we can reuse numerical
  # setups and we allocate one for each distinct aii.
  # * Otherwise, there will be no reuse so we only need one numerical setup
  # that will be updated.
  # To be general, we build a map sysslvrcaches: step -> NumericalSetup.
  # We will probably never need more than 256 stages so we can use Int8.
  if !reuse
    n = 1
    ptrs = fill(Int8(1), num_stages)
  else
    A = get_matrix(odeslvr.tableau)
    d = Dict{eltype(A),Int8}()
    n = 0
    ptrs = zeros(Int8, num_stages)
    for i in 1:num_stages
      aii = A[i, i]
      if !haskey(d, aii)
        n += 1
        d[aii] = n
      end
      ptrs[i] = d[aii]
    end
  end
  values = Vector{NumericalSetup}(undef, n)

  sysslvrcaches = CompressedArray(values, ptrs)
  odeslvrcache = (reuse, ui_pre, slopes, J, r, sysslvrcaches)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::DIMRungeKutta, odeop::ODEOperator{<:AbstractLinearODE},
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  reuse, ui_pre, slopes, J, r, sysslvrcaches = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr_l
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
    ws = (A[i, i] * dt, 1)

    # Update ODE operator cache
    update_odeopcache!(odeopcache, odeop, tx)

    # Create and solve stage operator
    # sysslvrcaches[i] will be unassigned at the first iteration
    sysslvrcache = isassigned(sysslvrcaches, i) ? sysslvrcaches[i] : nothing
    stageop = LinearStageOperator(
      odeop, odeopcache,
      tx, usx, ws,
      J, r, reuse, sysslvrcache
    )

    sysslvrcache = solve!(x, sysslvr, stageop, sysslvrcache)
    sysslvrcaches = _setindex_all!(sysslvrcaches, sysslvrcache, i)
  end

  # Update state
  tF = t0 + dt
  stateF = _update_dimrk!(stateF, state0, dt, slopes, b)

  # Pack outputs
  odeslvrcache = (reuse, ui_pre, slopes, J, r, sysslvrcaches)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _update_dimrk!(
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
