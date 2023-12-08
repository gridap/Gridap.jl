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
    The parameter ρ∞ of the generalized-α scheme must lie between zero and one.
    Setting ρ∞ to $(ρ∞01).
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
  t0::Real, us0::NTuple{2,AbstractVector}
)
  zero.(us0)
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha1,
  odeop::ODEOperator{<:AbstractLinearODE}, odeopcache,
  t0::Real, us0::NTuple{2,AbstractVector}
)
  J = allocate_jacobian(odeop, t0, us0, odeopcache)
  r = allocate_residual(odeop, t0, us0, odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  odeopcache, disopcache,
  u0::AbstractVector, v0::AbstractVector, dt::Real,
  tα::Real, αm::Real, αf::Real, γ::Real
)
  uα, vα = disopcache
  GeneralizedAlpha1NonlinearOperator(
    odeop, odeopcache,
    u0, v0, dt, tα, uα, vα, αm, αf, γ
  )
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator{<:AbstractLinearODE},
  odeopcache, disopcache,
  u0::AbstractVector, v0::AbstractVector, dt::Real,
  tα::Real, αm::Real, αf::Real, γ::Real
)
  J, r = disopcache
  GeneralizedAlpha1LinearOperator(
    odeop, odeopcache,
    u0, v0, dt, tα, αm, αf, γ,
    J, r
  )
end

function solve_odeop!(
  usF::NTuple{2,AbstractVector},
  odeslvr::GeneralizedAlpha1, odeop::ODEOperator,
  t0::Real, us0::NTuple{2,AbstractVector},
  cache
)
  u0, v0 = us0[1], us0[2]
  dt = get_dt(odeslvr)
  αm, αf, γ = odeslvr.αm, odeslvr.αf, odeslvr.γ
  tα = t0 + αf * dt
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
    u0, v0, dt, tα, αm, αf, γ
  )

  # Solve the discrete ODE operator
  usF, disslvrcache = solve_disop!(usF, odeslvr.disslvr, disop, disslvrcache)
  cache = (odeopcache, disopcache, disslvrcache)

  (tF, usF, cache)
end

######################
# Nonlinear operator #
######################
"""
    struct GeneralizedAlpha1NonlinearOperator <: DiscreteODEOperator end

Nonlinear discrete operator corresponding to the first-order generalized-α
scheme:
```math
residual(tx, ux, vx) = 0,

tx = (1 - αf) * t_n + αf * t_(n+1)
ux = (1 - αf) * u_n + αf * u_(n+1)
vx = (1 - αm) * v_n + αm * v_(n+1),

u_(n+1) = u_n + dt * ((1 - γ) * v_n + γ * x)
v_(n+1) = x.
```
"""
struct GeneralizedAlpha1NonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  u0::AbstractVector
  v0::AbstractVector
  dt::Real
  tα::Real
  uα::AbstractVector
  vα::AbstractVector
  αm::Real
  αf::Real
  γ::Real
end

function Algebra.residual!(
  r::AbstractVector,
  disop::GeneralizedAlpha1NonlinearOperator,
  x::AbstractVector
)
  u0, v0, dt = disop.u0, disop.v0, disop.dt
  tα, uα, vα = disop.tα, disop.uα, disop.vα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  # Residual: (uα, vα)
  uα, vα = _stage_alpha1!(uα, vα, x, u0, v0, dt, αm, αf, γ)

  tx = tα
  usx = (uα, vα)
  residual!(r, disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::GeneralizedAlpha1NonlinearOperator,
  x::AbstractVector
)
  tα, uα, vα = disop.tα, disop.uα, disop.vα
  tx = tα
  usx = (uα, vα)
  allocate_jacobian(disop.odeop, tx, usx, disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::GeneralizedAlpha1NonlinearOperator,
  x::AbstractVector
)
  u0, v0, dt = disop.u0, disop.v0, disop.dt
  tα, uα, vα = disop.tα, disop.uα, disop.vα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  # Jacobian: (uα, vα)
  uα, vα = _stage_alpha1!(uα, vα, x, u0, v0, dt, αm, αf, γ)

  tx = tα
  usx = (uα, vα)
  ws = (αf * γ * dt, αm)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tx, usx, ws, disop.odeopcache)
end

function solve_disop!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha1NonlinearOperator,
  disslvrcache
)
  uF, vF = usF[1], usF[2]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  u0, v0, dt = disop.u0, disop.v0, disop.dt
  tα = disop.tα
  γ = disop.γ

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, tα)

  # Solve the discrete ODE operator
  x = vF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalize
  uF, vF = _finalize_alpha1!(uF, vF, u0, v0, dt, γ)
  usF = (uF, vF)

  (usF, disslvrcache)
end

###################
# Linear operator #
###################
"""
    struct GeneralizedAlpha1LinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the first-order generalized-α scheme:
```math
residual(tx, ux, vx) = mass(tx) vx + stiffness(tx) ux + res(tx) = 0,

tx = (1 - αf) * t_n + αf * t_(n+1)
ux = (1 - αf) * u_n + αf * u_(n+1)
vx = (1 - αm) * v_n + αm * v_(n+1),

u_(n+1) = u_n + dt * ((1 - γ) * v_n + γ * x)
v_(n+1) = x.
```
"""
struct GeneralizedAlpha1LinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  u0::AbstractVector
  v0::AbstractVector
  dt::Real
  tα::Real
  αm::Real
  αf::Real
  γ::Real
  J::AbstractMatrix
  r::AbstractVector
end

# The jacobian matrix is constant if the ODE operator is linear and has
# constant forms
function is_jacobian_constant(disop::GeneralizedAlpha1LinearOperator)
  odeop = disop.odeop
  constant_jacobian = false
  if ODEOperatorType(odeop) <: AbstractLinearODE
    constant_stiffness = is_form_constant(odeop, 0)
    constant_mass = is_form_constant(odeop, 1)
    constant_jacobian = constant_stiffness && constant_mass
  end
  constant_jacobian
end

Algebra.get_matrix(disop::GeneralizedAlpha1LinearOperator) = disop.J

Algebra.get_vector(disop::GeneralizedAlpha1LinearOperator) = disop.r

function solve_disop!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlpha1LinearOperator,
  disslvrcache
)
  uF, vF = usF[1], usF[2]
  odeop, odeopcache = disop.odeop, disop.odeopcache
  u0, v0, dt = disop.u0, disop.v0, disop.dt
  tα = disop.tα
  αm, αf, γ = disop.αm, disop.αf, disop.γ
  J, r = disop.J, disop.r

  # Update the cache of the ODE operator (typically Dirichlet BCs)
  update_odeopcache!(odeopcache, odeop, tα)

  # Residual: (uα, vα)
  # Jacobian: (uα, vα)
  # Take x = 0 to split the mass term from the residual
  # Trick: use uF, vF to store uα, vα and r to store x = 0
  fill!(r, zero(eltype(r)))
  uα, vα = _stage_alpha1!(uF, vF, r, u0, v0, dt, αm, αf, γ)

  tx = tα
  usx = (uα, vα)
  ws = (αf * γ * dt, αm)
  # If the jacobian is constant, the following call will only retrieve the
  # jacobian from the ODEOpFromFEOpCache.
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, tx, usx, ws, disop.odeopcache)
  residual!(r, odeop, tx, usx, odeopcache)
  rmul!(r, -1)

  # Solve the discrete ODE operator
  x = vF
  disslvrcache = solve!(x, disslvr, disop, disslvrcache)

  # Finalize
  usF = _finalize_alpha1!(uF, vF, u0, v0, dt, γ)

  (usF, disslvrcache)
end

#########
# Utils #
#########
function _stage_alpha1!(
  ux::AbstractVector, vx::AbstractVector, x::AbstractVector,
  u0::AbstractVector, v0::AbstractVector, dt::Real,
  αm::Real, αf::Real, γ::Real
)
  # @. ux = (1 - αf) * u0 + αf * (u0 + dt * ((1 - γ) * v0 + γ * x))
  copy!(ux, u0)
  axpy!(αf * (1 - γ) * dt, v0, ux)
  axpy!(αf * γ * dt, x, ux)

  # @. vx = (1 - αm) * v0 + αm * x
  copy!(vx, v0)
  rmul!(vx, 1 - αm)
  axpy!(αm, x, vx)

  usx = (ux, vx)
  usx
end

function _finalize_alpha1!(
  uF::AbstractVector, vF::AbstractVector,
  u0::AbstractVector, v0::AbstractVector,
  dt::Real, γ::Real
)
  # @. uF = u0 + dt * ((1 - γ) * v0 + γ * x)
  copy!(uF, u0)
  axpy!((1 - γ) * dt, v0, uF)
  axpy!(γ * dt, vF, uF)

  usF = (uF, vF)
  usF
end
