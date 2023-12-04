"""
    struct IMEXRungeKutta <: ODESolver

Implicit-Explicit Runge-Kutta ODE solver.
"""
struct IMEXRungeKutta <: ODESolver
  disslvr_nl::NonlinearSolver
  disslvr_l::NonlinearSolver
  dt::Real
  tableau::AbstractTableau{ImplicitExplicitTableau}
end

# ODESolver interface
function get_dt(odeslvr::IMEXRungeKutta)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::ImplicitExplicitTableau,
  odeop::IMEXODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  im_ui, ex_ui = zero(x), zero(x)
  us = (x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (im_ui, ex_ui, J, r)
end

function allocate_disopcache(
  odeslvr::ImplicitExplicitTableau,
  odeop::IMEXODEOperator{LinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  us = (x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  im_res::AbstractVector{<:AbstractVector},
  ex_res::AbstractVector{<:AbstractVector},
  tableau::AbstractTableau
)
  im_ui, ex_ui, J, r = disopcache
  IMEXRungeKuttaNonlinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    im_ui, ex_ui, im_res, ex_res, tableau, t0, zero(t0), J, r
  )
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  im_res::AbstractVector{<:AbstractVector},
  ex_res::AbstractVector{<:AbstractVector},
  tableau::AbstractTableau
)
  J, r = disopcache
  IMEXRungeKuttaLinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    im_res, ex_res, tableau, J, r
  )
end

function solve_step!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, odeop::IMEXODEOperator,
  us0::NTuple{1,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, = us0
  dt = get_dt(odeslvr)
  tableau = get_tableau(odeslvr)
  im_tableau, ex_tableau = get_imex_tableaus(tableau)
  num_stages = length(get_nodes(im_tableau))

  # Allocate or unpack cache
  if isnothing(cache)
    us = (u0, u0)
    odeopcache = allocate_odeopcache(odeop, t0, us)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, u0)
    disslvrcache = allocate_disslvrcache(odeslvr)
    im_res = [zero(u0) for _ in 1:num_stages]
    ex_res = [zero(u0) for _ in 1:num_stages]
  else
    odeopcache, disopcache, disslvrcache, im_res, ex_res = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop,
    odeopcache, disopcache,
    t0, us0, dt,
    im_res, ex_res, tableau
  )

  # Solve the discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache, im_res, ex_res)

  (usF, tF, cache)
end

######################
# Nonlinear operator #
######################
"""
    struct IMEXRungeKuttaNonlinearOperator <: DiscreteODEOperator end

Nonlinear operator corresponding to a Implicit-Explicit Runge-Kutta scheme:
```math
residual(ti, ui, vi) = im_mass(ti, ui) vi
                     + ∑_{i < j} im_A[i, j] * im_res(tj, uj)
                     +           im_A[i, i] * im_res(ti, ui)
                     + ∑_{i < j} ex_A[i, j] * ex_res(tj, uj) = 0,

ti = t_n + c[i] * dt
ui = x
vi = (x - u_n) / dt,

residual(t_(n+1), u_(n+1), v_(n+1)) = im_mass(t_(n+1), u_(n+1)) v_(n+1)
                                    + ∑_{1 ≤ i ≤ s} im_b[i] * im_res(ti, ui)
                                    + ∑_{1 ≤ i ≤ s} ex_b[i] * ex_res(ti, ui) = 0,

u_(n+1) = x
v_(n+1) = (x - u_n) / dt.
```
"""
mutable struct IMEXRungeKuttaNonlinearOperator <: DiscreteODEOperator
  odeop::IMEXODEOperator
  odeopcache
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real

  im_res::AbstractVector{<:AbstractVector}
  ex_res::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  ti::Real
  aii::Real
  J::AbstractMatrix
  r::AbstractVector
end

# function Algebra.allocate_residual(
#   disop::IMEXRungeKuttaNonlinearOperator,
#   x::AbstractVector
# )
#   ti, dt, ui = disop.ti, disop.dt, disop.ui
#   !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
#   usi = (ui, x)
#   r = allocate_residual(disop.odeop, ti, usi, disop.odeopcache)
#   !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
#   r
# end

# function Algebra.residual!(
#   r::AbstractVector,
#   disop::SequentialRungeKuttaNonlinearOperator,
#   x::AbstractVector
# )
#   ti, dt, ui = disop.ti, disop.dt, disop.ui
#   !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
#   usi = (ui, x)
#   residual!(r, disop.odeop, ti, usi, disop.odeopcache)
#   !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
#   r
# end

# function Algebra.allocate_jacobian(
#   disop::SequentialRungeKuttaNonlinearOperator,
#   x::AbstractVector
# )
#   ti, dt, ui = disop.ti, disop.dt, disop.ui
#   !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
#   usi = (ui, x)
#   J = allocate_jacobian(disop.odeop, ti, usi, disop.odeopcache)
#   !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
#   J
# end

# function Algebra.jacobian!(
#   J::AbstractMatrix,
#   disop::SequentialRungeKuttaNonlinearOperator,
#   x::AbstractVector
# )
#   ti, dt, ui = disop.ti, disop.dt, disop.ui
#   ws = (disop.aii * dt, 1)
#   !iszero(disop.aii) && axpy!(disop.aii * dt, x, ui)
#   usi = (ui, x)
#   fillstored!(J, zero(eltype(J)))
#   jacobians!(J, disop.odeop, ti, usi, ws, disop.odeopcache)
#   !iszero(disop.aii) && axpy!(-disop.aii * dt, x, ui)
#   J
# end

# function Algebra.solve!(
#   usF::NTuple{1,AbstractVector},
#   odeslvr::RungeKutta, disop::SequentialRungeKuttaNonlinearOperator,
#   disslvrcaches
# )
# end

###################
# Linear operator #
###################
"""
    struct IMEXRungeKuttaLinearOperator <: LinearDiscreteODEOperator end

Linear operator corresponding to a Implicit-Explicit Runge-Kutta scheme:
```math
residual(ti, ui, vi) = im_mass(ti) vi
                     + ∑_{i < j} im_A[i, j] * (im_stiffness(tj) uj + im_res(tj))
                     +           im_A[i, i] * (im_stiffness(ti) ui + im_res(ti))
                     + ∑_{i < j} ex_A[i, j] * ex_res(tj, uj) = 0,

ti = t_n + c[i] * dt
ui = x
vi = (x - u_n) / dt,

residual(t_(n+1), u_(n+1), v_(n+1)) = im_mass(t_(n+1)) v_(n+1)
                                    + ∑_{1 ≤ i ≤ s} im_b[i] * (im_stiffness(ti) ui + im_res(ti))
                                    + ∑_{1 ≤ i ≤ s} ex_b[i] * ex_res(ti, ui) = 0,

u_(n+1) = x
v_(n+1) = (u_(n+1) - u_n) / dt.
```
"""
mutable struct IMEXRungeKuttaLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  us0::NTuple{1,AbstractVector}
  dt::Real
  im_res::AbstractVector{<:AbstractVector}
  ex_res::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::IMEXRungeKuttaLinearOperator) = disop.J

Algebra.get_vector(disop::IMEXRungeKuttaLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::RungeKutta, disop::IMEXRungeKuttaLinearOperator,
  disslvrcaches
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  t0, dt, us0 = disop.t0, disop.dt, disop.us0
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, res_temp = odeopcache
  im_res, ex_res = disop.im_res, disop.ex_res
  im_tableau, ex_tableau = get_imex_tableaus(disop.tableau)
  im_A, ex_A = get_matrix(im_tableau), get_matrix(ex_tableau)
  c = get_nodes(im_tableau)

  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  u0, = us0
  uF, = usF

  # Solve stages
  for i in eachindex(c)
    ti = t0 + c[i] * dt
    update_odeopcache!(odeopcache, odeop, ti)

    # Take linear combination of previous residuals
    fill!(res_temp, zero(eltype(res_temp)))
    for j in 1:i-1
      im_coef = im_A[i, j]
      if !iszero(im_coef)
        axpy!(im_coef, im_res[j], res_temp)
      end
      ex_coef = ex_A[i, j]
      if !iszero(ex_coef)
        axpy!(ex_coef, ex_res[j], res_temp)
      end
    end

    # Update jacobian and residual
    ui = uF
    fill!(ui, zero(eltype(ui)))
    vi = vs[i]
    copy!(vi, u0)
    rmul!(vi, -1 / dt)
    im_usi = (ui, vi) # ui does not matter
    ws = (im_A[i, i], 1 / dt)

    fillstored!(J, zero(eltype(J)))
    jacobians!(J, im_odeop, ti, im_usi, ws, im_odeopcache)

    filter = (true, false, false)
    residual!(r, im_odeop, ti, im_usi, im_odeopcache; filter)
    axpy!(im_A[i, i], r, res_temp)
    filter = (false, false, true)
    residual!(r, im_odeop, ti, im_usi, im_odeopcache; filter)
    axpy!(1, res_temp, r)
    rmul!(r, -1)

    # Solve the discrete ODE operator
    disslvrcache = solve!(ui, disslvr, disop, disslvrcache)

    # Store new residuals
    im_usi = (ui, ui) # vi does not matter
    ex_usi = (ui,)
    filter = (true, true, false)
    residual!(im_res[i], im_odeop, ti, usi, im_odeopcache; filter)
    residual!(ex_res[i], ex_odeop, ti, ex_usi, ex_odeopcache)
  end

  # Final mass system for the update
  usF, disslvrcache = _finalize_imex_rk!(usF, disslvr, disop, disslvrcache)

  disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)
  (usF, disslvrcaches)
end

#########
# Utils #
#########
function _finalize_imex_rk!(usF, disslvr, disop, disslvrcache)
  t0, dt, us0 = disop.t0, disop.dt, disop.us0
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  im_odeop, _ = get_imex_operators(odeop)
  im_odeopcache, _ = odeopcache

  im_res, ex_res = disop.im_res, disop.ex_res
  im_tableau, ex_tableau = get_imex_tableaus(disop.tableau)
  im_b = get_weights(im_tableau)
  ex_b = get_weights(ex_tableau)

  u0, = us0
  uF, = usF
  tF = t0 + dt

  im_usF = (uF, uF) # does not matter
  w = 1 / dt
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, im_odeop, tF, im_usF, 1, w, im_odeopcache)

  fill!(r, zero(eltype(r)))
  for i in eachindex(im_b, ex_b)
    im_coef = im_b[i]
    if !iszero(im_coef)
      axpy!(im_coef, im_res[i], r)
    end
    ex_coef = ex_b[i]
    if !iszero(ex_coef)
      axpy!(ex_coef, ex_res[i], r)
    end
  end
  axpy!(-1 / dt, u0, r)

  disslvrcache = solve!(uF, disslvr, disop, disslvrcache)
  (uF, disslvrcache)
end
