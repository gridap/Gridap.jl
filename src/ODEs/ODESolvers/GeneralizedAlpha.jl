"""
    struct GeneralizedAlpha <: ODESolver

Generalized-α ODE solver.
"""
struct GeneralizedAlpha <: ODESolver
  disslvr::NonlinearSolver
  dt::Float64
  ρ∞::Float64

  function GeneralizedAlpha(disslvr, dt, ρ∞)
    ρ∞01 = clamp(ρ∞, 0, 1)
    if ρ∞01 != ρ∞
      msg = """
      The parameter ρ∞ of the generalized-α must lie between zero and one.
      Setting ρ∞ to $(ρ∞01).
      """
      @warn msg
    end
    new(disslvr, dt, ρ∞01)
  end
end

# ODESolver interface
function get_dt(odeslvr::GeneralizedAlpha)
  odeslvr.dt
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha,
  odeop::ODEOperator, odeopcache,
  t::Real, x::AbstractVector
)
  ntuple(_ -> similar(x), get_order(odeop) + 1)
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha,
  odeop::ODEOperator{LinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  usx = ntuple(_ -> x, get_order(odeop) + 1)
  J = allocate_jacobian(odeop, t, usx, odeopcache)
  r = allocate_residual(odeop, t, usx, odeopcache)
  (J, r)
end

# 1st-order
function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tα::Real
)
  usα = disopcache
  GeneralizedAlphaNonlinearOperator1(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, usα
  )
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tα::Real
)
  J, r = disopcache
  GeneralizedAlphaLinearOperator1(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, J, r,
  )
end

function solve_step!(
  usF::NTuple{2,AbstractVector},
  odeslvr::GeneralizedAlpha, odeop::ODEOperator,
  us0::NTuple{2,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, v0 = us0
  dt = get_dt(odeslvr)
  ρ∞ = odeslvr.ρ∞

  αf = 1 / (1 + ρ∞)
  αm = (3 - ρ∞) / (1 + ρ∞) / 2
  γ = 1 / 2 + αm - αf
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

# 2nd-order
function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{3,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, β::Real, tα::Real
)
  usα = disopcache
  GeneralizedAlphaNonlinearOperator2(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, β, usα
  )
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{3,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, β::Real, tα::Real
)
  J, r = disopcache
  GeneralizedAlphaLinearOperator2(
    odeop, odeopcache,
    tα, us0, dt,
    αm, αf, γ, β, J, r
  )
end

function solve_step!(
  usF::NTuple{3,AbstractVector},
  odeslvr::GeneralizedAlpha, odeop::ODEOperator,
  us0::NTuple{3,AbstractVector}, t0::Real,
  cache
)
  # Unpack us and ODE solver
  u0, v0, a0 = us0
  dt = get_dt(odeslvr)
  ρ∞ = odeslvr.ρ∞

  αf = ρ∞ / (1 + ρ∞)
  αm = (2 * ρ∞ - 1) / (1 + ρ∞)
  γ = 1 / 2 - αm + αf
  β = (1 - αm + αf)^2 / 4
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

  # Solve discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (usF, tF, cache)
end

######################################
# GeneralizedAlphaNonlinearOperators #
######################################
"""
    struct GeneralizedAlphaNonlinearOperator1 <: DiscreteODEOperator end

Nonlinear discrete operator corresponding to the first-order generalized-α scheme:
```math
residual(tα, uα, vα) = 0,

uα = (1 - αf) * u_n + αf * u_(n+1)
vα = (1 - αm) * v_n + αm * v_(n+1),

u_(n+1) = u_n + h * ((1 - γ) * v_n + γ * x)
v_(n+1) = x.
```
"""
struct GeneralizedAlphaNonlinearOperator1 <: DiscreteODEOperator
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
  disop::GeneralizedAlphaNonlinearOperator1,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ)
  allocate_residual(disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::GeneralizedAlphaNonlinearOperator1,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ)
  residual!(r, disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::GeneralizedAlphaNonlinearOperator1,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ)
  allocate_jacobian(disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::GeneralizedAlphaNonlinearOperator1,
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
  disslvr::NonlinearSolver, disop::GeneralizedAlphaNonlinearOperator1,
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

"""
    struct GeneralizedAlphaNonlinearOperator2 <: DiscreteODEOperator end

Nonlinear discrete operator corresponding to the second-order generalized-α scheme:
```math
residual(tα, uα, vα, aα) = 0,

uα = αf * u_n + (1 - αf) * u_(n+1)
vα = αf * v_n + (1 - αf) * v_(n+1)
aα = αm * a_n + (1 - αm) * a_(n+1),

u_(n+1) = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
v_(n+1) = v0 + dt * ((1 - γ) * a0 + γ * x)
a_(n+1) = x.
```
"""
struct GeneralizedAlphaNonlinearOperator2 <: DiscreteODEOperator
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

function Algebra.allocate_residual(
  disop::GeneralizedAlphaNonlinearOperator2,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ, β)
  allocate_residual(disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::GeneralizedAlphaNonlinearOperator2,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ, β)
  residual!(r, disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::GeneralizedAlphaNonlinearOperator2,
  x::AbstractVector
)
  tα, dt, us0, usα = disop.tα, disop.dt, disop.us0, disop.usα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  usα = _fill_usα!(usα, us0, x, dt, αm, αf, γ, β)
  allocate_jacobian(disop.odeop, tα, usα, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::GeneralizedAlphaNonlinearOperator2,
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
  disslvr::NonlinearSolver, disop::GeneralizedAlphaNonlinearOperator2,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  tα, dt, us0 = disop.tα, disop.dt, disop.us0
  γ, β = disop.γ, disop.β

  update_odeopcache!(odeopcache, odeop, tα)

  uF, vF, aF = usF
  disslvrcache = solve!(aF, disslvr, disop, disslvrcache)
  usF = _fill_usF!(usF, us0, aF, dt, γ, β)

  (usF, disslvrcache)
end

###################################
# GeneralizedAlphaLinearOperators #
###################################
"""
    struct GeneralizedAlphaLinearOperator1 <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the first-order generalized-α scheme:
```math
residual(tα, uα, vα) = mass(tα) vα + stiffness(tαf) uα + res(tαf) = 0,

uα = (1 - αf) * u_n + αf * u_(n+1)
vα = (1 - αm) * v_n + αm * v_(n+1),

u_(n+1) = u_n + h * ((1 - γ) * v_n + γ * x)
v_(n+1) = x.
```
"""
struct GeneralizedAlphaLinearOperator1 <: LinearDiscreteODEOperator
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

Algebra.get_matrix(disop::GeneralizedAlphaLinearOperator1) = disop.J

Algebra.get_vector(disop::GeneralizedAlphaLinearOperator1) = disop.r

function Algebra.solve!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlphaLinearOperator1,
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

"""
    struct GeneralizedAlphaLinearOperator2 <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the second-order generalized-α scheme:
```math
residual(tα, uα, vα, aα) = mass(tα) aα + damping(tα) vα + stiffness(tα) uα + res(tα) = 0,

uα = αf * u_n + (1 - αf) * u_(n+1)
vα = αf * v_n + (1 - αf) * v_(n+1)
aα = αm * a_n + (1 - αm) * a_(n+1),

u_(n+1) = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
v_(n+1) = v0 + dt * ((1 - γ) * a0 + γ * x)
a_(n+1) = x.
```
"""
struct GeneralizedAlphaLinearOperator2 <: LinearDiscreteODEOperator
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

Algebra.get_matrix(disop::GeneralizedAlphaLinearOperator2) = disop.J

Algebra.get_vector(disop::GeneralizedAlphaLinearOperator2) = disop.r

function Algebra.solve!(
  usF::NTuple{3,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlphaLinearOperator2,
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

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tα, usα, ws, disop.odeopcache)
  residual!(r, odeop, tα, usα, odeopcache)
  rmul!(r, -1)

  aF = usF[3]
  disslvrcache = solve!(aF, disslvr, disop, disslvrcache)
  usF = _fill_usF!(usF, us0, aF, dt, γ, β)

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

function _fill_usF!(
  usF::NTuple{3,AbstractVector}, us0::NTuple{3,AbstractVector}, x,
  dt, γ, β
)
  u0, v0, a0 = us0
  uF, vF, aF = usF
  @. uF = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
  @. vF = v0 + dt * ((1 - γ) * a0 + γ * x)
  @. aF = x
  (uF, vF, aF)
end

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
