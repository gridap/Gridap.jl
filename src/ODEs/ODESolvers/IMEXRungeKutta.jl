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
  t0::Real, us0::NTuple{2,AbstractVector}
)
  throw(imex_rk_not_implemented_msg)
end

function allocate_disopcache(
  odeslvr::IMEXRungeKutta,
  odeop::IMEXODEOperator{<:AbstractSemilinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  u0 = us0[1]
  ui, ri = zero(u0), zero(u0)
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (ui, ri, J, r)
end

function allocate_disopcache(
  odeslvr::IMEXRungeKutta,
  odeop::IMEXODEOperator{<:AbstractLinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (J, r)
end

function allocate_disslvrcache(odeslvr::IMEXRungeKutta)
  (nothing, nothing)
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator,
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
  im_res::AbstractVector{<:AbstractVector},
  ex_res::AbstractVector{<:AbstractVector},
  tableau::AbstractTableau
)
  throw(imex_rk_not_implemented_msg)
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator{<:AbstractSemilinearODE},
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
  im_res::AbstractVector{<:AbstractVector},
  ex_res::AbstractVector{<:AbstractVector},
  tableau::AbstractTableau
)
  ui, ri, J, r = disopcache
  ti = t0
  aii = zero(t0)
  IMEXRungeKuttaNonlinearOperator(
    odeop, odeopcache,
    t0, u0, dt,
    ti, ui, aii, ri,
    im_res, ex_res, tableau,
    J, r
  )
end

function DiscreteODEOperator(
  odeslvr::IMEXRungeKutta, odeop::IMEXODEOperator{<:AbstractLinearODE},
  odeopcache, disopcache,
  t0::Real, u0::AbstractVector, dt::Real,
  im_res::AbstractVector{<:AbstractVector},
  ex_res::AbstractVector{<:AbstractVector},
  tableau::AbstractTableau
)
  J, r = disopcache
  IMEXRungeKuttaLinearOperator(
    odeop, odeopcache,
    t0, u0, dt,
    im_res, ex_res, tableau,
    J, r
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
  t0::Real, us0::NTuple{1,AbstractVector},
  cache
)
  u0 = us0[1]
  dt = get_dt(odeslvr)
  tF = t0 + dt

  tableau = odeslvr.tableau
  im_tableau, _ = get_imex_tableaus(tableau)
  num_stages = length(get_nodes(im_tableau))

  # Allocate or unpack cache
  if isnothing(cache)
    us0_full = (u0, u0)
    odeopcache = allocate_odeopcache(odeop, t0, us0_full)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, us0_full)
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
    t0, u0, dt,
    im_res, ex_res, tableau
  )

  # Solve the discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr, disop, disslvrcache)
  cache = (odeopcache, disopcache, disslvrcache, im_res, ex_res)

  (tF, usF, cache)
end

######################
# Nonlinear operator #
######################
"""
    struct IMEXRungeKuttaNonlinearOperator <: DiscreteODEOperator end

Nonlinear operator corresponding to an Implicit-Explicit Runge-Kutta scheme:
```math
residual(tx, ux, vx) = im_mass(tx) vx
                     + ∑_{i < j} im_A[i, j] * im_res[j]
                     +           im_A[i, i] * im_res(tx, ux)
                     + ∑_{i < j} ex_A[i, j] * ex_res[j] = 0,

tx = t_n + c[i] * dt
ux = u_n + dt * x
vx = x,

residual(t_(n+1), u_(n+1), v_(n+1)) = im_mass(t_(n+1)) v_(n+1)
                                    + ∑_{1 ≤ i ≤ s} im_b[i] * im_res[i]
                                    + ∑_{1 ≤ i ≤ s} ex_b[i] * ex_res[i] = 0,

u_(n+1) = u_n + dt * x
v_(n+1) = x.
```
"""
mutable struct IMEXRungeKuttaNonlinearOperator <: DiscreteODEOperator
  odeop::IMEXODEOperator
  odeopcache
  t0::Real
  u0::AbstractVector
  dt::Real
  ti::Real
  ui::AbstractVector
  aii::Real
  ri::AbstractVector
  im_res::AbstractVector{<:AbstractVector}
  ex_res::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::AbstractMatrix
  r::AbstractVector
end

function Algebra.residual!(
  r::AbstractVector,
  disop::IMEXRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  u0, dt = disop.u0, disop.dt
  ti, ui, aii, ri = disop.ti, disop.ui, disop.aii, disop.ri

  im_odeop, _ = get_imex_operators(disop.odeop)
  im_odeopcache, _, tmp = disop.odeopcache

  # Pre weight the massless residual
  # Residual: (u0 + dt * x, 0)
  # Trick: use tmp to store 0
  copy!(ui, u0)
  axpy!(dt, x, ui)
  fill!(tmp, zero(eltype(tmp)))

  tx = ti
  usx = (ui, tmp)
  residual!(r, im_odeop, tx, usx, im_odeopcache)
  rmul!(r, aii - 1)

  # Full residual
  # Residual: (u0 + dt * x, x)
  tx = ti
  usx = (ui, x)
  residual!(tmp, im_odeop, tx, usx, im_odeopcache)
  axpy!(1, tmp, r)

  # Residual: linear combination of previous residuals
  axpy!(1, ri, r)

  r
end

function Algebra.allocate_jacobian(
  disop::IMEXRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  ti, ui = disop.ti, disop.ui

  im_odeop, _ = get_imex_operators(disop.odeop)
  im_odeopcache, _ = disop.odeopcache

  tx = ti
  usx = (ui, x)
  allocate_jacobian(im_odeop, tx, usx, im_odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::IMEXRungeKuttaNonlinearOperator,
  x::AbstractVector
)
  u0, dt = disop.u0, disop.dt
  ti, ui, aii = disop.ti, disop.ui, disop.aii

  im_odeop, _ = get_imex_operators(disop.odeop)
  im_odeopcache, _ = disop.odeopcache

  # Jacobian: (u0 + dt * x, x)
  copy!(ui, u0)
  axpy!(dt, x, ui)

  tx = ti
  usx = (ui, x)
  ws = (aii * dt, 1)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, im_odeop, tx, usx, ws, im_odeopcache)
end

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, disop::IMEXRungeKuttaNonlinearOperator,
  disslvrcaches
)
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, u0, dt = disop.t0, disop.u0, disop.dt

  ri = disop.ri
  im_res, ex_res = disop.im_res, disop.ex_res
  tableau = disop.tableau
  im_tableau, ex_tableau = get_imex_tableaus(disop.tableau)
  im_A, ex_A = get_matrix(im_tableau), get_matrix(ex_tableau)
  c = get_nodes(im_tableau)

  J, r = disop.J, disop.r

  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, _ = odeopcache

  # Solve stages
  for i in eachindex(c)
    ti = t0 + c[i] * dt

    # Update the cache of the ODE operator (typically Dirichlet BCs)
    update_odeopcache!(odeopcache, odeop, ti)

    # Update operator state
    disop.ti = ti
    disop.aii = im_A[i, i]

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

    # Solve the discrete ODE operator
    explicit = iszero(disop.aii)
    disslvr, islvr = get_solver_index(odeslvr, explicit)
    disslvrcache = disslvrcaches[islvr]

    # Decide whether the operator is linear or nonlinear
    if explicit
      # Jacobian: (any, any)
      tx = ti
      usx = (u0, u0)
      w = 1
      fillstored!(J, zero(eltype(J)))
      jacobian!(J, im_odeop, tx, usx, 1, w, im_odeopcache)

      # Residual: linear combination of the stored residuals
      copy!(r, ri)
      rmul!(r, -1)

      _op = IMEXRungeKuttaLinearOperator(
        odeop, odeopcache,
        t0, u0, dt,
        im_res, ex_res, tableau, J, r
      )
    else
      _op = disop
    end

    x = uF
    disslvrcache = solve!(x, disslvr, _op, disslvrcache)
    disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)

    # Store residuals at current stage
    # Residual: (u0 + dt * x, 0)
    # Trick: use r to store 0
    axpby!(1, u0, dt, x)
    ux = x
    fill!(r, zero(eltype(r)))

    tx = ti
    im_usx = (ux, r)
    ex_usx = (ux,)
    residual!(im_res[i], im_odeop, tx, im_usx, im_odeopcache)
    residual!(ex_res[i], ex_odeop, tx, ex_usx, ex_odeopcache)
  end

  # Finalize
  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  finalop = IMEXRungeKuttaFinalizingOperator(
    odeop, odeopcache,
    t0, u0, dt,
    im_res, ex_res, tableau, J, r
  )
  usF, disslvrcaches = solve!(usF, odeslvr, finalop, disslvrcaches)

  (usF, disslvrcaches)
end

###################
# Linear operator #
###################
"""
    struct IMEXRungeKuttaLinearOperator <: LinearDiscreteODEOperator end

Linear operator corresponding to an Implicit-Explicit Runge-Kutta scheme:
```math
residual(tx, ux, vx) = im_mass(tx) vx
                     + ∑_{i < j} im_A[i, j] * im_res[j]
                     +           im_A[i, i] * (im_stiffness(tx) ux + im_res(tx))
                     + ∑_{i < j} ex_A[i, j] * ex_res[j] = 0,

tx = t_n + c[i] * dt
ux = u_n + dt * x
vx = x,

residual(t_(n+1), u_(n+1), v_(n+1)) = im_mass(t_(n+1)) v_(n+1)
                                    + ∑_{1 ≤ i ≤ s} im_b[i] * im_res[i]
                                    + ∑_{1 ≤ i ≤ s} ex_b[i] * ex_res[i] = 0,

u_(n+1) = un + dt * x
v_(n+1) = x.
```
"""
mutable struct IMEXRungeKuttaLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  u0::AbstractVector
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
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, u0, dt = disop.t0, disop.u0, disop.dt

  im_res, ex_res = disop.im_res, disop.ex_res
  tableau = disop.tableau
  im_tableau, ex_tableau = get_imex_tableaus(tableau)
  im_A, ex_A = get_matrix(im_tableau), get_matrix(ex_tableau)
  c = get_nodes(im_tableau)

  J, r = disop.J, disop.r

  im_odeop, ex_odeop = get_imex_operators(odeop)
  im_odeopcache, ex_odeopcache, _ = odeopcache

  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  # Solve stages
  for i in eachindex(c)
    ti = t0 + c[i] * dt

    # Update the cache of the ODE operator (typically Dirichlet BCs)
    update_odeopcache!(odeopcache, odeop, ti)

    # Jacobian: (any, any)
    tx = ti
    usx = (u0, u0)
    ws = (im_A[i, i] * dt, 1)
    fillstored!(J, zero(eltype(J)))
    jacobians!(J, im_odeop, tx, usx, ws, im_odeopcache)

    # Implicit residual
    # Residual: (u0, 0)
    # Trick: use uF to store 0
    fill!(uF, zero(eltype(uF)))

    tx = ti
    usx = (u0, uF)
    residual!(r, im_odeop, tx, usx, im_odeopcache)
    rmul!(r, im_A[i, i])

    # Residual: linear combination of residuals at previous stages
    for j in 1:i-1
      im_coef = im_A[i, j]
      if !iszero(im_coef)
        axpy!(im_coef, im_res[j], r)
      end
      ex_coef = ex_A[i, j]
      if !iszero(ex_coef)
        axpy!(ex_coef, ex_res[j], r)
      end
    end
    rmul!(r, -1)

    # Solve the discrete ODE operator
    x = uF
    disslvrcache = solve!(x, disslvr, disop, disslvrcache)

    # Store residuals at current stage
    # (u0 + dt * x, 0)
    # Trick: use r to store 0
    axpby!(1, u0, dt, x)
    ux = x
    fill!(r, zero(eltype(r)))
    im_usx = (ux, r)
    ex_usx = (ux,)
    residual!(im_res[i], im_odeop, tx, im_usx, im_odeopcache)
    residual!(ex_res[i], ex_odeop, tx, ex_usx, ex_odeopcache)
  end

  # Finalizing step
  finalop = IMEXRungeKuttaFinalizingOperator(
    odeop, odeopcache,
    t0, u0, dt,
    im_res, ex_res, tableau, J, r
  )
  usF, disslvrcaches = solve!(usF, odeslvr, finalop, disslvrcaches)

  (usF, disslvrcaches)
end

#############
# Finalizer #
#############
"""
    struct IMEXRungeKuttaFinalizingOperator <: LinearDiscreteODEOperator end

Linear operator corresponding to the final update of an Implicit-Explicit
Runge-Kutta scheme:
```math
residual(t_(n+1), u_(n+1), v_(n+1)) = im_mass(t_(n+1)) v_(n+1)
                                    + ∑_{1 ≤ i ≤ s} im_b[i] * im_res[i]
                                    + ∑_{1 ≤ i ≤ s} ex_b[i] * ex_res[i] = 0,

u_(n+1) = u_n + dt * x
v_(n+1) = x.
```
"""
mutable struct IMEXRungeKuttaFinalizingOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  t0::Real
  u0::AbstractVector
  dt::Real
  im_res::AbstractVector{<:AbstractVector}
  ex_res::AbstractVector{<:AbstractVector}
  tableau::AbstractTableau
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::IMEXRungeKuttaFinalizingOperator) = disop.J

Algebra.get_vector(disop::IMEXRungeKuttaFinalizingOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{1,AbstractVector},
  odeslvr::IMEXRungeKutta, disop::IMEXRungeKuttaFinalizingOperator,
  disslvrcaches
)
  uF = usF[1]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, u0, dt = disop.t0, disop.u0, disop.dt

  im_res, ex_res = disop.im_res, disop.ex_res
  im_tableau, ex_tableau = get_imex_tableaus(disop.tableau)
  im_b, ex_b = get_weights(im_tableau), get_weights(ex_tableau)

  J, r = disop.J, disop.r

  im_odeop, _ = get_imex_operators(odeop)
  im_odeopcache, _ = odeopcache

  # Jacobian: (any, any)
  tx = t0 + dt
  usx = (u0, u0)
  w = 1
  fillstored!(J, zero(eltype(J)))
  jacobian!(J, im_odeop, tx, usx, 1, w, im_odeopcache)

  # Residual: linear combination of previous residuals
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
  rmul!(r, -1)

  explicit = true
  disslvr, islvr = get_solver_index(odeslvr, explicit)
  disslvrcache = disslvrcaches[islvr]

  x = uF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)
  disslvrcaches = Base.setindex(disslvrcaches, disslvrcache, islvr)

  # @. uF = u0 + dt * x
  axpby!(1, u0, dt, x)
  uF = x
  usF = (uF,)

  (usF, disslvrcaches)
end
