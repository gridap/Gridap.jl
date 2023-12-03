"""
    struct GeneralizedAlpha1 <: ODESolver

Generalized-α first-order ODE solver.
"""
struct GeneralizedAlpha1 <: ODESolver
  disslvr::NonlinearSolver
  dt::Float64
  αm::Float64
  αf::Float64
  γ::Float64
end

# Constructors
function GeneralizedAlpha1(
  disslvr::NonlinearSolver,
  dt::Float64, ρ∞::Float64
)
  ρ∞01 = clamp(ρ∞, 0, 1)
  if ρ∞01 != ρ∞
    msg = """
    The parameter ρ∞ of the generalized-α scheme must lie between zero and one. Setting ρ∞ to $(ρ∞01).
    """
    @warn msg
    ρ∞ = ρ∞01
  end

  αm = (3 - ρ∞) / (1 + ρ∞) / 2
  αf = 1 / (1 + ρ∞)
  γ = 1 / 2 + αm - αf
  GeneralizedAlpha1(disslvr, dt, αm, αf, γ)
end

# ODESolver interface
function get_dt(odeslvr::GeneralizedAlpha1)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha1,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  (similar(x), similar(x),)
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha1,
  odeop::ODEOperator{LinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  usx = (x, x)
  J = allocate_jacobian(odeop, t, usx, odeopcache)
  r = allocate_residual(odeop, t, usx, odeopcache)
  (J, r)
end

# 1st-order
function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tα::Real
)
  usα = disopcache
  GeneralizedAlpha1NonlinearOperator(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, usα
  )
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tα::Real
)
  J, r = disopcache
  GeneralizedAlpha1LinearOperator(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, J, r,
  )
end

function solve_step!(
  usF::NTuple{2,AbstractVector},
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  us0::NTuple{2,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, v0 = us0
  dt = get_dt(odeslvr)
  αm, αf, γ = odeslvr.αm, odeslvr.αf, odeslvr.γ
  tα = t0 + αf * dt

  # Allocate or unpack cache
  if isnothing(cache)
    odeopcache = allocate_odeopcache(odeop, t0, us0)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, u0)
    disslvrcache = allocate_disslvrcache(odeslvr)
  else
    odeopcache, disopcache, disslvrcache = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop,
    odeopcache, disopcache,
    t0, us0, dt,
    αm, αf, γ, tα
  )

  # Solve discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (usF, tF, cache)
end

######################
# Nonlinear operator #
######################
"""
    struct GeneralizedAlpha1NonlinearOperator <: DiscreteODEOperator end

Nonlinear discrete operator corresponding to the first-order generalized-α scheme:
```math
residual(tα, uα, vα) = 0,

uα = (1 - αf) * u_n + αf * u_(n+1)
vα = (1 - αm) * v_n + αm * v_(n+1),

u_(n+1) = u_n + h * ((1 - γ) * v_n + γ * x)
v_(n+1) = x.
```
"""
struct GeneralizedAlpha1NonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  tα::Real
  us0::NTuple{2,AbstractVector}
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  usα::NTuple{2,AbstractVector}
end

function Algebra.allocate_residual(
  disop::GeneralizedAlpha1NonlinearOperator,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ)
  allocate_residual(disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::GeneralizedAlpha1NonlinearOperator,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ)
  residual!(r, disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::GeneralizedAlpha1NonlinearOperator,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ)
  allocate_jacobian(disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::GeneralizedAlpha1NonlinearOperator,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ)
  ws = _get_ws(dt, αm, αf, γ)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tα, usα, ws, disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha1NonlinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  tα, dt, us0 = disop.tα, disop.dt, disop.us0
  γ = disop.γ

  update_odeopcache!(odeopcache, odeop, tα)

  uF, vF = usF
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  usF = _fill_usF!(usF, us0, vF, dt, γ)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct GeneralizedAlpha1LinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the first-order generalized-α scheme:
```math
residual(tα, uα, vα) = mass(tα) vα + stiffness(tαf) uα + res(tαf) = 0,

uα = (1 - αf) * u_n + αf * u_(n+1)
vα = (1 - αm) * v_n + αm * v_(n+1),

u_(n+1) = u_n + h * ((1 - γ) * v_n + γ * x)
v_(n+1) = x.
```
"""
struct GeneralizedAlpha1LinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  tα::Real
  us0::NTuple{2,AbstractVector}
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::GeneralizedAlpha1LinearOperator) = disop.J

Algebra.get_vector(disop::GeneralizedAlpha1LinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha1LinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  tα, dt, us0 = disop.tα, disop.dt, disop.us0
  αf, αm, γ = disop.αf, disop.αm, disop.γ

  update_odeopcache!(odeopcache, odeop, tα)

  u0, v0 = us0
  uα, vα = usF
  @. uα = (1 - αf) * u0 + αf * (u0 + dt * (1 - γ) * v0)
  @. vα = (1 - αm) * v0
  usα = (uα, vα)
  ws = _get_ws(dt, αm, αf, γ)

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tα, usα, ws, disop.odeopcache)
  residual!(r, odeop, tα, usα, odeopcache)
  rmul!(r, -1)

  vF = usF[2]
  disslvrcache = solve!(vF, disslvr, disop, disslvrcache)
  usF = _fill_usF!(usF, us0, vF, dt, γ)

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _fill_usF!(
  usF::NTuple{2,AbstractVector}, us0::NTuple{2,AbstractVector}, x,
  dt, γ
)
  u0, v0 = us0
  uF, vF = usF
  @. uF = u0 + dt * ((1 - γ) * v0 + γ * x)
  @. vF = x
  (uF, vF)
end

function _get_ws(dt, αm, αf, γ)
  wu = αf * dt * γ
  wv = αm
  (wu, wv)
end

function _fill_usα!(
  usα::NTuple{2,AbstractVector}, us0::NTuple{2,AbstractVector}, x,
  dt, αm, αf, γ
)
  u0, v0 = us0
  uα, vα = usα
  @. uα = (1 - αf) * u0 + αf * (u0 + dt * ((1 - γ) * v0 + γ * x))
  @. vα = (1 - αm) * v0 + αm * x
  (uα, vα)
end
