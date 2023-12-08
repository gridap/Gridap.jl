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
  t0::Real, us0::NTuple{3,AbstractVector}
)
  zero.(us0)
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha2,
  odeop::ODEOperator{<:AbstractLinearODE}, odeopcache,
  t0::Real, us0::NTuple{3,AbstractVector}
)
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  odeopcache, disopcache,
  u0::AbstractVector, v0::AbstractVector, a0::AbstractVector, dt::Real,
  tα::Real, αm::Real, αf::Real, γ::Real, β::Real
)
  uα, vα, aα = disopcache
  GeneralizedAlpha2NonlinearOperator(
    odeop, odeopcache,
    u0, v0, a0, dt,
    tα, uα, vα, aα, αm, αf, γ, β
  )
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator{<:AbstractLinearODE},
  odeopcache, disopcache,
  u0::AbstractVector, v0::AbstractVector, a0::AbstractVector, dt::Real,
  tα::Real, αm::Real, αf::Real, γ::Real, β::Real,
)
  J, r = disopcache
  GeneralizedAlpha2LinearOperator(
    odeop, odeopcache,
    u0, v0, a0, dt,
    tα, αm, αf, γ, β, J, r
  )
end

function solve_odeop!(
  usF::NTuple{3,AbstractVector},
  odeslvr::GeneralizedAlpha2, odeop::ODEOperator,
  t0::Real, us0::NTuple{3,AbstractVector},
  cache
)
  u0, v0, a0 = us0
  dt = get_dt(odeslvr)
  αm, αf, γ, β = odeslvr.αm, odeslvr.αf, odeslvr.γ, odeslvr.β
  tα = t0 + (1 - αf) * dt
  tF = t0 + dt

  # Allocate or unpack cache
  if isnothing(cache)
    odeopcache = allocate_odeopcache(odeop, t0, us0)
    disopcache = allocate_disopcache(odeslvr, odeop, odeopcache, t0, us0)
    disslvrcache = allocate_disslvrcache(odeslvr)
  else
    odeopcache, disopcache, disslvrcache = cache
  end

  # Create discrete ODE operator
  disop = DiscreteODEOperator(
    odeslvr, odeop,
    odeopcache, disopcache,
    u0, v0, a0, dt,
    tα, αm, αf, γ, β
  )

  # Solve the discrete ODE operator
  usF, disslvrcache = solve_disop!(usF, odeslvr.disslvr, disop, disslvrcache)

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (tF, usF, cache)
end

######################
# Nonlinear operator #
######################
"""
    struct GeneralizedAlpha2NonlinearOperator <: DiscreteODEOperator end

Nonlinear discrete operator corresponding to the second-order generalized-α scheme:
```math
residual(tx, ux, vx, ax) = 0,

tx = αf * t_n + (1 - αf) * t_(n+1)
ux = αf * u_n + (1 - αf) * u_(n+1)
vx = αf * v_n + (1 - αf) * v_(n+1)
ax = αm * a_n + (1 - αm) * a_(n+1),

u_(n+1) = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
v_(n+1) = v0 + dt * ((1 - γ) * a0 + γ * x)
a_(n+1) = x.
```
"""
struct GeneralizedAlpha2NonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  u0::AbstractVector
  v0::AbstractVector
  a0::AbstractVector
  dt::Real
  tα::Real
  uα::AbstractVector
  vα::AbstractVector
  aα::AbstractVector
  αm::Real
  αf::Real
  γ::Real
  β::Real
end

function Algebra.residual!(
  r::AbstractVector,
  disop::GeneralizedAlpha2NonlinearOperator,
  x::AbstractVector
)
  u0, v0, a0, dt = disop.u0, disop.v0, disop.a0, disop.dt
  tα, uα, vα, aα = disop.tα, disop.uα, disop.vα, disop.aα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  # Residual: (uα, vα, aα)
  uα, vα, aα = _stage_alpha2!(uα, vα, aα, x, u0, v0, a0, dt, αm, αf, γ, β)

  tx = tα
  usx = (uα, vα, aα)
  residual!(r, disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::GeneralizedAlpha2NonlinearOperator,
  x::AbstractVector
)
  tα, uα, vα, aα = disop.tα, disop.uα, disop.vα, disop.aα
  tx = tα
  usx = (uα, vα, aα)
  allocate_jacobian(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::GeneralizedAlpha2NonlinearOperator,
  x::AbstractVector
)
  u0, v0, a0, dt = disop.u0, disop.v0, disop.a0, disop.dt
  tα, uα, vα, aα = disop.tα, disop.uα, disop.vα, disop.aα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  # Jacobian: (uα, vα, aα)
  uα, vα, aα = _stage_alpha2!(uα, vα, aα, x, u0, v0, a0, dt, αm, αf, γ, β)

  tx = tα
  usx = (uα, vα, aα)
  ws = ((1 - αf) * β * dt^2, (1 - αf) * γ * dt, 1 - αm)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tx, usx, ws, disop.odeopcache)
end

function solve_disop!(
  usF::NTuple{3,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha2NonlinearOperator,
  disslvrcache
)
  uF, vF, aF = usF[1], usF[2], usF[3]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  u0, v0, a0, dt = disop.u0, disop.v0, disop.a0, disop.dt
  tα = disop.tα
  γ, β = disop.γ, disop.β

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, tα)

  # Solve the discrete ODE operator
  x = aF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalize
  uF, vF, aF = _finalize_alpha2!(uF, vF, aF, u0, v0, a0, dt, γ, β)
  usF = (uF, vF, aF)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct GeneralizedAlpha2LinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the second-order generalized-α scheme:
```math
residual(tx, ux, vx, ax) = mass(tx) ax + damping(tx) vx + stiffness(tx) ux + res(tx) = 0,

tx = αf * t_n + (1 - αf) * t_(n+1)
ux = αf * u_n + (1 - αf) * u_(n+1)
vx = αf * v_n + (1 - αf) * v_(n+1)
ax = αm * a_n + (1 - αm) * a_(n+1),

u_(n+1) = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
v_(n+1) = v0 + dt * ((1 - γ) * a0 + γ * x)
a_(n+1) = x.
```
"""
struct GeneralizedAlpha2LinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  u0::AbstractVector
  v0::AbstractVector
  a0::AbstractVector
  dt::Real
  tα::Real
  αm::Real
  αf::Real
  γ::Real
  β::Real
  J::AbstractMatrix
  r::AbstractVector
end

# The jacobian matrix is constant if the ODE operator is linear and has
# constant forms
function is_jacobian_constant(disop::GeneralizedAlpha2LinearOperator)
  odeop = disop.odeop
  constant_jacobian = false
  if ODEOperatorType(odeop) <: AbstractLinearODE
    constant_stiffness = is_form_constant(odeop, 0)
    constant_damping = is_form_constant(odeop, 1)
    constant_mass = is_form_constant(odeop, 2)
    constant_jacobian = constant_stiffness && constant_damping && constant_mass
  end
  constant_jacobian
end

Algebra.get_matrix(disop::GeneralizedAlpha2LinearOperator) = disop.J

Algebra.get_vector(disop::GeneralizedAlpha2LinearOperator) = disop.r

function solve_disop!(
  usF::NTuple{3,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha2LinearOperator,
  disslvrcache
)
  uF, vF, aF = usF[1], usF[2], usF[3]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  u0, v0, a0, dt = disop.u0, disop.v0, disop.a0, disop.dt
  tα = disop.tα
  αm, αf, γ, β = disop.αm, disop.αf, disop.γ, disop.β
  J, r = disop.J, disop.r

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, tα)

  # Residual: (uα, vα, aα)
  # Jacobian: (uα, vα, aα)
  # Take x = 0 to split the mass term from the residual
  # Trick: use uF, vF, aF to store uα, vα, aα and r to store x = 0
  fill!(r, zero(eltype(r)))
  uα, vα, aα = _stage_alpha2!(uF, vF, aF, r, u0, v0, a0, dt, αm, αf, γ, β)

  tx = tα
  usx = (uα, vα, aα)
  ws = ((1 - αf) * β * dt^2, (1 - αf) * γ * dt, 1 - αm)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tx, usx, ws, disop.odeopcache)
  residual!(r, odeop, tx, usx, odeopcache)
  rmul!(r, -1)

  # Solve the discrete ODE operator
  x = aF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Express usF in terms of the solution of the discrete ODE operator
  uF, vF, aF = _finalize_alpha2!(uF, vF, aF, u0, v0, a0, dt, γ, β)
  usF = (uF, vF, aF)

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _stage_alpha2!(
  ux::AbstractVector, vx::AbstractVector, ax::AbstractVector, x::AbstractVector,
  u0::AbstractVector, v0::AbstractVector, a0::AbstractVector, dt::Real,
  αm::Real, αf::Real, γ::Real, β::Real
)
  # @. ux = αf * u0 + (1 - αf) * (u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x))
  copy!(ux, u0)
  axpy!((1 - αf) * dt, v0, ux)
  axpy!((1 - αf) * (1 - 2 * β) * dt^2 / 2, a0, ux)
  axpy!((1 - αf) * β * dt^2, x, ux)

  # @. vx = αf * v0 + (1 - αf) * (v0 + dt * ((1 - γ) * a0 + γ * x))
  copy!(vx, v0)
  axpy!((1 - αf) * (1 - γ) * dt, a0, vx)
  axpy!((1 - αf) * γ * dt, x, vx)

  # @. ax = αm * a0 + (1 - αm) * x
  copy!(ax, a0)
  rmul!(ax, αm)
  axpy!(1 - αm, x, ax)

  usx = (ux, vx, ax)
  usx
end

function _finalize_alpha2!(
  uF::AbstractVector, vF::AbstractVector, aF::AbstractVector,
  u0::AbstractVector, v0::AbstractVector, a0::AbstractVector, dt::Real,
  γ::Real, β::Real
)
  # @. uF = u0 + dt * v0 + dt^2 / 2 * ((1 - 2 * β) * a0 + 2 * β * x)
  copy!(uF, u0)
  axpy!(dt, v0, uF)
  axpy!((1 - 2 * β) * dt^2 / 2, a0, uF)
  axpy!(β * dt^2, aF, uF)

  # @. vF = v0 + dt * ((1 - γ) * a0 + γ * x)
  copy!(vF, v0)
  axpy!((1 - γ) * dt, a0, vF)
  axpy!(γ * dt, aF, vF)

  usF = (uF, vF, aF)
  usF
end
