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

const imex_rk_not_implemented_msg = """
IMEX Runge-Kutta is only implemented for IMEX ODE operators whose implicit
residual is semilinear.
"""

function allocate_disopcache(
  odeslvr::IMEXRungeKutta,
  odeop::IMEXODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  throw(imex_rk_not_implemented_msg)
end

function allocate_disopcache(
  odeslvr::IMEXRungeKutta,
  odeop::IMEXODEOperator{<:AbstractSemilinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  vi, ri = zero(x), zero(x)
  us = (x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (vi, ri, J, r)
end

function allocate_disopcache(
  odeslvr::IMEXRungeKutta,
  odeop::IMEXODEOperator{<:AbstractLinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  us = (x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (J, r)
end

function allocate_disslvrcache(odeslvr::IMEXRungeKutta)
  (nothing, nothing)
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  im_res::AbstractVector{<:AbstractVector},
  ex_res::AbstractVector{<:AbstractVector},
  tableau::AbstractTableau
)
  throw(imex_rk_not_implemented_msg)
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator{<:AbstractSemilinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
  im_res::AbstractVector{<:AbstractVector},
  ex_res::AbstractVector{<:AbstractVector},
  tableau::AbstractTableau
)
  vi, ri, J, r = disopcache
  ti, aii = t0, zero(t0)
  IMEXRungeKuttaNonlinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    im_res, ex_res, tableau,
    ti, aii, vi, ri, J, r
  )
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator{<:AbstractLinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{1,AbstractVector}, dt::Real,
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

# RungeKutta interface
function get_tableau(odeslvr::IMEXRungeKutta)
  odeslvr.tableau
end

function get_solver_index(odeslvr::IMEXRungeKutta, explicit::Bool)
  explicit ? (odeslvr.disslvr_l, 2) : (odeslvr.disslvr_nl, 1)
end

function solve_step!(
  usF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator,
  us0::NTuple{1,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, = us0
  dt = get_dt(odeslvr)
  tableau = odeslvr.tableau
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
residual(ti, ui, vi) = im_mass(ti) vi
                     + ∑_{i < j} im_A[i, j] * im_res(tj, uj)
                     +           im_A[i, i] * im_res(ti, ui)
                     + ∑_{i < j} ex_A[i, j] * ex_res(tj, uj) = 0,

ti = t_n + c[i] * dt
ui = x
vi = (x - u_n) / dt,

residual(t_(n+1), u_(n+1), v_(n+1)) = im_mass(t_(n+1)) v_(n+1)
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
  vi::AbstractVector
  ri::AbstractVector
  J::AbstractMatrix
  r::AbstractVector
end

function Algebra.residual!(
  r::AbstractVector,
  disop::IMEXRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, dt, (u0,) = disop.ti, disop.dt, disop.us0
  vi, ri, aii = disop.vi, disop.ri, disop.aii
  im_odeop, _ = get_imex_operators(disop.odeop)
  im_odeopcache, _, res_temp = disop.odeopcache

  @. vi = (x - u0) / dt
  usi = (x, vi)
  # Mass matrix
  filter = (false, false, true)
  residual!(r, im_odeop, ti, usi, im_odeopcache; filter)
  # Implicit res
  filter = (true, true, false)
  residual!(res_temp, im_odeop, ti, usi, im_odeopcache; filter)
  axpy!(aii, res_temp, r)
  # Precomputed linear combination of previous residuals
  axpy!(1, ri, r)
  r
end

function Algebra.allocate_jacobian(
  disop::IMEXRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti = disop.ti
  im_odeop, _ = get_imex_operators(disop.odeop)
  im_odeopcache, _ = disop.odeopcache

  usi = (x, x) # vi does not matter
  allocate_jacobian(im_odeop, ti, usi, im_odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::IMEXRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, dt = disop.ti, disop.dt
  aii = disop.aii
  im_odeop, _ = get_imex_operators(disop.odeop)
  im_odeopcache, _ = disop.odeopcache

  usi = (x, x) # vi does not matter
  ws = (aii, 1 / dt)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, im_odeop, ti, usi, ws, im_odeopcache)
  J
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, disop::IMEXRungeKuttaNonlinearOperator,
  disslvrcaches
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, _ = odeopcache

  t0, dt, us0 = disop.t0, disop.dt, disop.us0
  J, r = disop.J, disop.r, ri = disop.J, disop.r, disop.ri

  im_res, ex_res = disop.im_res, disop.ex_res
  tableau = disop.tableau
  im_tableau, ex_tableau = get_imex_tableaus(disop.tableau)
  im_A, ex_A = get_matrix(im_tableau), get_matrix(ex_tableau)
  c = get_nodes(im_tableau)

  u0, = us0
  uF, = usF

  # Solve stages
  for i in eachindex(c)
    ti = t0 + c[i] * dt
    update_odeopcache!(odeopcache, odeop, ti)

    # Take linear combination of previous residuals
    fill!(ri, zero(eltype(ri)))
    for j in 1:i-1
      im_coef = im_A[i, j]
      if !iszero(im_coef)
        axpy!(im_coef, im_res[j], ri)
      end
      ex_coef = ex_A[i, j]
      if !iszero(ex_coef)
        axpy!(ex_coef, ex_res[j], ri)
      end
    end

    # Update operator state
    disop.ti = ti
    disop.aii = im_A[i, i]

    # Solve the discrete ODE operator
    explicit = iszero(im_A[i, i])
    disslvr, islvr = get_solver_index(odeslvr, explicit)
    disslvrcache = disslvrcaches[islvr]

    if explicit
      copy!(uF, u0)
      rmul!(uF, -1 / dt)
      usi = (uF, uF) # ui does not matter
      filter = (false, false, true)
      w = 1 / dt

      fillstored!(J, zero(eltype(J)))
      jacobian!(J, im_odeop, ti, usi, 1, w, im_odeopcache)
      residual!(r, im_odeop, ti, usi, im_odeopcache; filter)
      axpy!(1, ri, r)
      rmul!(r, -1)

      _op = IMEXRungeKuttaLinearOperator(
        odeop, odeopcache,
        t0, us0, dt,
        im_res, ex_res, tableau, J, r
      )
    else
      _op = disop
    end

    ui = uF
    disslvrcache = solve!(ui, disslvr, _op, disslvrcache)
    disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)

    # Store new residuals
    im_usi = (ui, ui)
    ex_usi = (ui,)
    filter = (true, true, false)
    residual!(im_res[i], im_odeop, ti, im_usi, im_odeopcache; filter)
    residual!(ex_res[i], ex_odeop, ti, ex_usi, ex_odeopcache; filter)
  end

  # Final mass system for the update
  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  _op = IMEXRungeKuttaLinearOperator(
    odeop, odeopcache,
    t0, us0, dt,
    im_res, ex_res, tableau, J, r
  )
  usF, disslvrcache = _finalize_imex_rk!(usF, disslvr, _op, disslvrcache)

  disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)
  (usF, disslvrcaches)
end

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
  odeslvr::IMEXRungeKutta, disop::IMEXRungeKuttaLinearOperator,
  disslvrcaches
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, res_temp = odeopcache

  t0, dt, us0 = disop.t0, disop.dt, disop.us0
  J, r = disop.J, disop.r

  im_res, ex_res = disop.im_res, disop.ex_res
  tableau = disop.tableau
  im_tableau, ex_tableau = get_imex_tableaus(tableau)
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
    ui, vi = uF, uF
    copy!(vi, u0)
    rmul!(vi, -1 / dt)
    im_usi = (ui, vi) # ui does not matter
    ws = (im_A[i, i], 1 / dt)

    fillstored!(J, zero(eltype(J)))
    jacobians!(J, im_odeop, ti, im_usi, ws, im_odeopcache)
    # Extract rest of residual
    filter = (true, false, false)
    residual!(r, im_odeop, ti, im_usi, im_odeopcache; filter)
    axpy!(im_A[i, i], r, res_temp)
    # Extract rest of mass matrix
    filter = (false, false, true)
    residual!(r, im_odeop, ti, im_usi, im_odeopcache; filter)
    axpy!(1, res_temp, r)
    rmul!(r, -1)

    # Solve the discrete ODE operator
    disslvrcache = solve!(ui, disslvr, disop, disslvrcache)

    # Store new residuals
    im_usi = (ui, ui)
    ex_usi = (ui,)
    filter = (true, true, false)
    residual!(im_res[i], im_odeop, ti, im_usi, im_odeopcache; filter)
    residual!(ex_res[i], ex_odeop, ti, ex_usi, ex_odeopcache; filter)
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
  J, r = disop.J, disop.r
  odeop, odeopcache = disop.odeop, disop.odeopcache
  im_odeop, _ = get_imex_operators(odeop)
  im_odeopcache, _ = odeopcache

  im_res, ex_res = disop.im_res, disop.ex_res
  im_tableau, ex_tableau = get_imex_tableaus(disop.tableau)
  im_b, ex_b = get_weights(im_tableau), get_weights(ex_tableau)

  u0, = us0
  uF, = usF
  tF = t0 + dt

  im_usF = (uF, uF) # uF does not matter but vF does
  copy!(uF, u0)
  rmul!(uF, -1 / dt)
  w = 1 / dt
  filter = (false, false, true)

  fillstored!(J, zero(eltype(J)))
  jacobian!(J, im_odeop, tF, im_usF, 1, w, im_odeopcache)
  residual!(r, im_odeop, tF, im_usF, im_odeopcache; filter)

  # Linear combination of residuals
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
  rmul!(r, -1)

  disslvrcache = solve!(uF, disslvr, disop, disslvrcache)

  usF = (uF,)
  (usF, disslvrcache)
end
