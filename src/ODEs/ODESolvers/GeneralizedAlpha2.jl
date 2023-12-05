"""
    struct GeneralizedAlpha2 <: ODESolver

Generalized-α second-order ODE solver.
"""
struct GeneralizedAlpha2 <: ODESolver
  disslvr::NonlinearSolver
  dt::Float64
  αm::Float64
  αf::Float64
  γ::Float64
  β::Float64
end

# Constructors
function GeneralizedAlpha2(disslvr::NonlinearSolver, dt::Float64, ρ∞::Float64)
  ρ∞01 = clamp(ρ∞, 0, 1)
  if ρ∞01 != ρ∞
    msg = """
    The parameter ρ∞ of the generalized-α scheme must lie between zero and one.
    Setting ρ∞ to $(ρ∞01).
    """
    @warn msg
    ρ∞ = ρ∞01
  end

  αf = ρ∞ / (1 + ρ∞)
  αm = (2 * ρ∞ - 1) / (1 + ρ∞)
  γ = 1 / 2 - αm + αf
  β = (1 - αm + αf)^2 / 4
  GeneralizedAlpha2(disslvr, dt, αm, αf, γ, β)
end

function Newmark(disslvr::NonlinearSolver, dt::Float64, γ::Float64, β::Float64)
  γ01 = clamp(γ, 0, 1)
  if γ01 != γ
    msg = """
    The parameter γ of the Newmark scheme must lie between zero and one.
    Setting γ to $(γ01).
    """
    @warn msg
    γ = γ01
  end

  β01 = clamp(β, 0, 1)
  if β01 != β
    msg = """
    The parameter β of the Newmark scheme must lie between zero and one.
    Setting β to $(β01).
    """
    @warn msg
    β = β01
  end

  αm, αf = 0.0, 0.01
  GeneralizedAlpha2(disslvr, dt, αm, αf, γ, β)
end

# ODESolver interface
function get_dt(odeslvr::GeneralizedAlpha2)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha2,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  (zero(x), zero(x), zero(x))
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha2,
  odeop::ODEOperator{<:AbstractLinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  us = (x, x, x)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  r = allocate_residual(odeop, t, us, odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{3,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, β::Real, tα::Real
)
  usα = disopcache
  GeneralizedAlpha2NonlinearOperator(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, β, usα
  )
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator{<:AbstractLinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{3,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, β::Real, tα::Real
)
  J, r = disopcache
  GeneralizedAlpha2LinearOperator(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, β, J, r
  )
end

function solve_step!(
  usF::NTuple{3,AbstractVector},
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  us0::NTuple{3,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, v0, a0 = us0
  dt = get_dt(odeslvr)
  αm, αf, γ, β = odeslvr.αm, odeslvr.αf, odeslvr.γ, odeslvr.β
  tα = t0 + (1 - αf) * dt

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
    αm, αf, γ, β, tα
  )

  # Solve the discrete ODE operator
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
    struct GeneralizedAlpha2NonlinearOperator <: DiscreteODEOperator end

Nonlinear discrete operator corresponding to the second-order generalized-α scheme:
```math
residual(tα, uα, vα, aα) = 0,

tα = αf * t_n + (1 - αf) * t_(n+1)
uα = αf * u_n + (1 - αf) * u_(n+1)
vα = αf * v_n + (1 - αf) * v_(n+1)
aα = αm * a_n + (1 - αm) * a_(n+1),

u_(n+1) = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
v_(n+1) = v0 + dt * ((1 - γ) * a0 + γ * x)
a_(n+1) = x.
```
"""
struct GeneralizedAlpha2NonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  tα::Real
  us0::NTuple{3,AbstractVector}
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  β::Real
  usα::NTuple{3,AbstractVector}
end

function Algebra.residual!(
  r::AbstractVector,
  disop::GeneralizedAlpha2NonlinearOperator,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ, β)
  residual!(r, disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::GeneralizedAlpha2NonlinearOperator,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ, β)
  allocate_jacobian(disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::GeneralizedAlpha2NonlinearOperator,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ, β)
  ws = _get_ws(dt, αm, αf, γ, β)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tα, usα, ws, disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{3,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha2NonlinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  tα, dt, us0 = disop.tα, disop.dt, disop.us0
  γ, β = disop.γ, disop.β

  update_odeopcache!(odeopcache, odeop, tα)

  uF, vF, aF = usF
  disslvrcache = solve!(aF, disslvr, disop, disslvrcache)

  # Express usF in terms of the solution of the discrete ODE operator
  usF = _finalize_alpha2!(usF, us0, aF, dt, γ, β)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct GeneralizedAlpha2LinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the second-order generalized-α scheme:
```math
residual(tα, uα, vα, aα) = mass(tα) aα + damping(tα) vα + stiffness(tα) uα + res(tα) = 0,

tα = αf * t_n + (1 - αf) * t_(n+1)
uα = αf * u_n + (1 - αf) * u_(n+1)
vα = αf * v_n + (1 - αf) * v_(n+1)
aα = αm * a_n + (1 - αm) * a_(n+1),

u_(n+1) = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
v_(n+1) = v0 + dt * ((1 - γ) * a0 + γ * x)
a_(n+1) = x.
```
"""
struct GeneralizedAlpha2LinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  tα::Real
  us0::NTuple{3,AbstractVector}
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  β::Real
  J::AbstractMatrix
  r::AbstractVector
end

Algebra.get_matrix(disop::GeneralizedAlpha2LinearOperator) = disop.J

Algebra.get_vector(disop::GeneralizedAlpha2LinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{3,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha2LinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  J, r = disop.J, disop.r
  tα, dt, us0 = disop.tα, disop.dt, disop.us0
  αf, αm, γ, β = disop.αf, disop.αm, disop.γ, disop.β

  update_odeopcache!(odeopcache, odeop, tα)

  u0, v0, a0 = us0
  uα, vα, aα = usF
  @. aα = αm * a0
  @. vα = αf * v0 + (1 - αf) * (v0 + dt * (1 - γ) * a0)
  @. uα = αf * u0 + (1 - αf) * (u0 + dt * v0 + dt^2 / 2 * (1 - 2 * β) * a0)
  usα = (uα, vα, aα)
  ws = _get_ws(dt, αm, αf, γ, β)

  # Update jacobian and residual
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tα, usα, ws, disop.odeopcache)
  residual!(r, odeop, tα, usα, odeopcache)
  rmul!(r, -1)

  # Solve the discrete ODE operator
  aF = usF[3]
  disslvrcache = solve!(aF, disslvr, disop, disslvrcache)

  # Express usF in terms of the solution of the discrete ODE operator
  usF = _finalize_alpha2!(usF, us0, aF, dt, γ, β)

  (usF, disslvrcache)
end

#############
# Finalizer #
#############
function _finalize_alpha2!(
  usF::NTuple{3,AbstractVector}, us0::NTuple{3,AbstractVector},
  x::AbstractVector, dt::Real, γ::Real, β::Real
)
  u0, v0, a0 = us0
  uF, vF, aF = usF
  @. uF = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
  @. vF = v0 + dt * ((1 - γ) * a0 + γ * x)
  @. aF = x
  (uF, vF, aF)
end

#########
# Utils #
#########
function _get_ws(dt, αm, αf, γ, β)
  wu = (1 - αf) * dt^2 * β
  wv = (1 - αf) * dt * γ
  wa = 1 - αm
  (wu, wv, wa)
end

function _fill_usα!(
  usα::NTuple{3,AbstractVector}, us0::NTuple{3,AbstractVector}, x,
  dt, αm, αf, γ, β
)
  u0, v0, a0 = us0
  uα, vα, aα = usα
  @. uα = αf * u0 + (1 - αf) * (u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x))
  @. vα = αf * v0 + (1 - αf) * (v0 + dt * ((1 - γ) * a0 + γ * x))
  @. aα = αm * a0 + (1 - αm) * x
  (uα, vα, aα)
end
