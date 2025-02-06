"""
    struct IMEXRungeKutta <: ODESolver

Implicit-Explicit Runge-Kutta ODE solver.
```
mass(tx, ux) vx + im_res(tx, ux) = 0,

tx = t_n + c[i] * dt
ux = u_n + ∑_{1 ≤ j < i} im_A[i, j] * dt * im_slopes[j] + im_A[i, i] * dt * x
         + ∑_{1 ≤ j < i} ex_A[i, j] * dt * ex_slopes[j]
vx = x
im_slopes[i] = x,

mass(tx, ux) vx + ex_res(tx, ux) = 0,

tx = t_n + c[i] * dt
ux = u_n + ∑_{1 ≤ j ≤ i} im_A[i, j] * dt * im_slopes[j]
         + ∑_{1 ≤ j < i} ex_A[i, j] * dt * ex_slopes[j]
vx = x
ex_slopes[i] = x,

u_(n+1) = u_n + ∑_{1 ≤ i ≤ s} im_b[i] * dt * im_slopes[i]
              + ∑_{1 ≤ i ≤ s} ex_b[i] * dt * ex_slopes[i].
```
"""
struct IMEXRungeKutta <: ODESolver
  sysslvr_nl::NonlinearSolver
  sysslvr_l::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{ImplicitExplicitTableau}
end

#######################
# Notimplemented case #
#######################
const imex_rk_not_implemented_msg = """
IMEX Runge-Kutta is only implemented for IMEX ODE operators whose implicit
residual is quasilinear.
"""

function allocate_odecache(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  @unreachable imex_rk_not_implemented_msg
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  @unreachable imex_rk_not_implemented_msg
end

# Dispatch on the IMEX decomposition
function allocate_odecache(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator{<:AbstractQuasilinearODE},
  t0::Real, us0::NTuple{1,AbstractVector}
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  allocate_odecache(odeslvr, odeop, im_odeop, ex_odeop, t0, us0)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator{<:AbstractQuasilinearODE},
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  im_odeop, ex_odeop = get_imex_operators(odeop)
  ode_march!(stateF, odeslvr, odeop, im_odeop, ex_odeop, t0, state0, odecache)
end

##################
# Nonlinear case #
##################
# This is very similar to `DIMRungeKutta` applied to a nonlinear `ODEOperator`
function allocate_odecache(
  odeslvr::IMEXRungeKutta, odeop,
  im_odeop::ODEOperator{<:AbstractQuasilinearODE},
  ex_odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  ui_pre, ui = zero(u0), zero(u0)
  im_tableau, ex_tableau = get_imex_tableaus(odeslvr.tableau)
  im_num_stages = length(get_nodes(im_tableau))
  im_slopes = [zero(u0) for _ in 1:im_num_stages]
  ex_num_stages = length(get_nodes(ex_tableau))
  ex_slopes = [zero(u0) for _ in 1:ex_num_stages]

  is_semilinear = ODEOperatorType(im_odeop) <: AbstractSemilinearODE
  mass_constant = is_form_constant(odeop, 1)
  reuse = (is_semilinear && mass_constant)

  # The explicit stage is always going to be linear so we always need to
  # allocate a jacobian and residual.
  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  # We can share the numerical setup across the implicit and explicit parts
  # because the linear solver will only be called when the implicit part goes
  # through an explicit stage (aii = 0) and on the explicit part. In both cases
  # the matrix of the linear stage operator is the mass matrix.
  sysslvrcaches = (nothing, nothing)
  odeslvrcache = (reuse, ui_pre, ui, im_slopes, ex_slopes, J, r, sysslvrcaches)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, odeop,
  im_odeop::ODEOperator{<:AbstractQuasilinearODE},
  ex_odeop::ODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  im_odeopcache, ex_odeopcache = odeopcache
  reuse, ui_pre, ui, im_slopes, ex_slopes, J, r, sysslvrcaches = odeslvrcache
  sysslvrcache_nl, sysslvrcache_l = sysslvrcaches

  # Unpack solver
  sysslvr_nl, sysslvr_l = odeslvr.sysslvr_nl, odeslvr.sysslvr_l
  dt, tableau = odeslvr.dt, odeslvr.tableau
  im_tableau, ex_tableau = get_imex_tableaus(tableau)
  im_A, im_b = get_matrix(im_tableau), get_weights(im_tableau)
  ex_A, ex_b = get_matrix(ex_tableau), get_weights(ex_tableau)
  c = get_nodes(im_tableau)

  for i in eachindex(c)
    # Define scheme
    x = im_slopes[i]
    tx = t0 + c[i] * dt
    copy!(ui_pre, u0)
    for j in 1:i-1
      axpy!(im_A[i, j] * dt, im_slopes[j], ui_pre)
      axpy!(ex_A[i, j] * dt, ex_slopes[j], ui_pre)
    end

    # Update ODE operator cache
    update_odeopcache!(odeopcache, odeop, tx)

    # 1. Implicit part
    aii = im_A[i, i]

    # If the implicit tableau is padded, we can skip the first implicit solve
    # and set im_slopes[1] to zero
    if i == 1 && is_padded(tableau)
      fill!(x, zero(eltype(x)))
    else
      # The stage becomes explicit when aii = 0 because the implicit part is
      # quasilinear
      explicit_stage = iszero(aii)

      if explicit_stage
        # Define scheme
        # Set x to zero to split jacobian and residual
        fill!(x, zero(eltype(x)))
        usx = (ui_pre, x)
        ws = (0, 1)

        # Create and solve stage operator
        im_stageop = LinearStageOperator(
          im_odeop, im_odeopcache,
          tx, usx, ws,
          J, r, reuse, sysslvrcache_l
        )

        sysslvrcache_l = solve!(x, sysslvr_l, im_stageop, sysslvrcache_l)
      else
        # Define scheme
        function usx(x)
          copy!(ui, ui_pre)
          axpy!(aii * dt, x, ui)
          (ui, x)
        end
        ws = (aii * dt, 1)

        # Create and solve stage operator
        im_stageop = NonlinearStageOperator(
          im_odeop, im_odeopcache,
          tx, usx, ws
        )

        sysslvrcache_nl = solve!(x, sysslvr_nl, im_stageop, sysslvrcache_nl)
      end
    end

    # 2. Explicit part
    # vx does not matter
    # Compute ui from ui_pre and im_slopes[i]
    copy!(ui, ui_pre)
    axpy!(aii * dt, im_slopes[i], ui)
    usx = (ui, u0)

    # This stage operator is a little more complicated than usual so we build
    # the jacobian and residual here:
    # [m(ti, ui)] x + [ex_res(ti, ui)] = 0
    # * The explicit part does not contain the mass, we take the full residual
    # * The jacobian is the mass matrix, which is stored in the implicit part
    residual!(r, ex_odeop, tx, usx, ex_odeopcache)
    if isnothing(sysslvrcache_l) || !reuse
      ws = (0, 1)
      jacobian!(J, im_odeop, tx, usx, ws, im_odeopcache)
    end
    ex_stageop = LinearStageOperator(J, r, reuse)

    x = ex_slopes[i]
    sysslvrcache_l = solve!(x, sysslvr_l, ex_stageop, sysslvrcache_l)
  end

  # Update state
  tF = t0 + dt
  stateF = _update_imexrk!(stateF, state0, dt, im_slopes, im_b, ex_slopes, ex_b)

  # Pack outputs
  sysslvrcaches = (sysslvrcache_nl, sysslvrcache_l)
  odeslvrcache = (reuse, ui_pre, ui, im_slopes, ex_slopes, J, r, sysslvrcaches)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

###############
# Linear case #
###############
# This is very similar to `DIMRungeKutta` applied to a linear `ODEOperator`
function allocate_odecache(
  odeslvr::IMEXRungeKutta, odeop,
  im_odeop::ODEOperator{<:AbstractLinearODE},
  ex_odeop::ODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = allocate_odeopcache(odeop, t0, us0N)

  ui_pre = zero(u0)
  im_tableau, ex_tableau = get_imex_tableaus(odeslvr.tableau)
  im_num_stages = length(get_nodes(im_tableau))
  im_slopes = [zero(u0) for _ in 1:im_num_stages]
  ex_num_stages = length(get_nodes(ex_tableau))
  ex_slopes = [zero(u0) for _ in 1:ex_num_stages]

  stiffness_constant = is_form_constant(im_odeop, 0)
  mass_constant = is_form_constant(im_odeop, 1)
  reuse = (stiffness_constant && mass_constant)

  # The explict part will always bring about a linear stage operator, and its
  # matrix is always the mass matrix. From the constraint on the nodes of the
  # implicit and explicit tableaus, we know that the first stage of the
  # implicit part is in fact explicit so the corresponding matrix is also the
  # mass matrix. This means that the numerical setup of the explict part can
  # always be chosen as the same as numerical setup of the first stage of the
  # implicit part.

  # For the same reason, the jacobian and residual can be shared across the
  # implicit and explicit parts.
  J = allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = allocate_residual(odeop, t0, us0N, odeopcache)

  # Numerical setups for the linear solver
  # * If the mass and stiffness matrices are constant, we can reuse numerical
  # setups and we allocate one for each distinct aii.
  # * Otherwise, there will be no reuse so we only need one numerical setup
  # that will be updated.
  # To be general, we build a map sysslvrcaches: step -> NumericalSetup.
  # We will probably never need more than 256 stages so we can use Int8
  if !reuse
    n = 1
    ptrs = fill(Int8(1), im_num_stages)
  else
    im_A = get_matrix(im_tableau)
    d = Dict{eltype(im_A),Int8}()
    n = 0
    ptrs = zeros(Int8, im_num_stages)
    for i in 1:im_num_stages
      aii = im_A[i, i]
      if !haskey(d, aii)
        n += 1
        d[aii] = n
      end
      ptrs[i] = d[aii]
    end
  end
  values = Vector{NumericalSetup}(undef, n)

  sysslvrcaches = CompressedArray(values, ptrs)
  odeslvrcache = (reuse, ui_pre, im_slopes, ex_slopes, J, r, sysslvrcaches)

  (odeslvrcache, odeopcache)
end

function ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, odeop,
  im_odeop::ODEOperator{<:AbstractLinearODE},
  ex_odeop::ODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache = odecache
  im_odeopcache, ex_odeopcache = odeopcache
  reuse, ui_pre, im_slopes, ex_slopes, J, r, sysslvrcaches = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr_l
  dt, tableau = odeslvr.dt, odeslvr.tableau
  im_tableau, ex_tableau = get_imex_tableaus(tableau)
  im_A, im_b = get_matrix(im_tableau), get_weights(im_tableau)
  ex_A, ex_b = get_matrix(ex_tableau), get_weights(ex_tableau)
  c = get_nodes(im_tableau)

  for i in eachindex(c)
    # Define scheme
    # Set x to zero to split jacobian and residual
    x = im_slopes[i]
    fill!(x, zero(eltype(x)))
    tx = t0 + c[i] * dt
    copy!(ui_pre, u0)
    for j in 1:i-1
      axpy!(im_A[i, j] * dt, im_slopes[j], ui_pre)
      axpy!(ex_A[i, j] * dt, ex_slopes[j], ui_pre)
    end
    usx = (ui_pre, x)
    ws = (im_A[i, i] * dt, 1)

    # Update ODE operator cache
    update_odeopcache!(odeopcache, odeop, tx)

    # Create and solve stage operator
    # 1. Implicit stage
    # If the implicit tableau is padded, we can skip the first implicit solve
    # and set im_slopes[1] to zero
    if i == 1 && is_padded(tableau)
      fill!(x, zero(eltype(x)))
    else
      # sysslvrcaches[i] will be unassigned at the first iteration
      sysslvrcache = isassigned(sysslvrcaches, i) ? sysslvrcaches[i] : nothing

      im_stageop = LinearStageOperator(
        im_odeop, im_odeopcache,
        tx, usx, ws,
        J, r, reuse, sysslvrcache
      )

      sysslvrcache = solve!(x, sysslvr, im_stageop, sysslvrcache)
      sysslvrcaches = _setindex_all!(sysslvrcaches, sysslvrcache, i)
    end

    # 2. Explicit part
    # vx does not matter
    # Compute ui from ui_pre and im_slopes[i]
    ui = axpy!(im_A[i, i] * dt, im_slopes[i], ui_pre)
    usx = (ui, u0)

    # This stage operator is a little more complicated than usual so we build
    # the jacobian and residual here:
    # [m(ti, ui)] x + [ex_res(ti, ui)] = 0
    # * The explicit part does not contain the mass, we take the full residual
    # * The jacobian is the mass matrix, which is stored in the implicit part
    residual!(r, ex_odeop, tx, usx, ex_odeopcache)
    sysslvrcache = isassigned(sysslvrcaches, 1) ? sysslvrcaches[1] : nothing
    if isnothing(sysslvrcache) || !reuse
      ws = (0, 1)
      jacobian!(J, im_odeop, tx, usx, ws, im_odeopcache)
    end
    ex_stageop = LinearStageOperator(J, r, reuse)

    x = ex_slopes[i]
    sysslvrcache = solve!(x, sysslvr, ex_stageop, sysslvrcache)
    sysslvrcaches = _setindex_all!(sysslvrcaches, sysslvrcache, i)
  end

  # Update state
  tF = t0 + dt
  stateF = _update_imexrk!(stateF, state0, dt, im_slopes, im_b, ex_slopes, ex_b)

  # Pack outputs
  odeslvrcache = (reuse, ui_pre, im_slopes, ex_slopes, J, r, sysslvrcaches)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache)
end

#########
# Utils #
#########
function _update_imexrk!(
  stateF::NTuple{1,AbstractVector}, state0::NTuple{1,AbstractVector},
  dt::Real, im_slopes::AbstractVector{<:AbstractVector}, im_b::AbstractVector,
  ex_slopes::AbstractVector{<:AbstractVector}, ex_b::AbstractVector
)
  # uF = u0 + ∑_{1 ≤ i ≤ s} im_b[i] * dt * im_slopes[i]
  #         + ∑_{1 ≤ i ≤ s} ex_b[i] * dt * ex_slopes[i]
  u0 = state0[1]
  uF = stateF[1]

  copy!(uF, u0)
  for (bi, slopei) in zip(im_b, im_slopes)
    axpy!(bi * dt, slopei, uF)
  end
  for (bi, slopei) in zip(ex_b, ex_slopes)
    axpy!(bi * dt, slopei, uF)
  end

  (uF,)
end
