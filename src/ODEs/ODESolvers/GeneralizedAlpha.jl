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
  similar(x)
end

function allocate_disopcache(
  odeslvr::GeneralizedAlpha,
  odeop::ODEOperator{LinearODE}, odeopcache,
  t::Real, x::AbstractVector
)
  J = allocate_jacobian(odeop, t, (x, x), odeopcache)
  r = allocate_residual(odeop, t, (x, x), odeopcache)
  (J, r)
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha, odeop::ODEOperator,
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tαf::Real
)
  vαm = disopcache
  GeneralizedAlphaNonlinearOperator(
    odeop, odeopcache, vαm,
    t0, us0, dt, αm, αf, γ, tαf
  )
end

function DiscreteODEOperator(
  odeslvr::GeneralizedAlpha, odeop::ODEOperator{LinearODE},
  odeopcache, disopcache,
  t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tαf::Real
)
  J, r = disopcache
  GeneralizedAlphaLinearOperator(
    odeop, odeopcache, J, r,
    t0, us0, dt, αm, αf, γ, tαf
  )
end

###############
# solve_step! #
###############
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
  tαf = t0 + αf * dt

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
    αm, αf, γ, tαf
  )

  # Solve discrete ODE operator
  usF, disslvrcache = solve!(usF, odeslvr.disslvr, disop, disslvrcache)
  tF = t0 + dt

  # Update cache
  cache = (odeopcache, disopcache, disslvrcache)

  (usF, tF, cache)
end

#####################################
# GeneralizedAlphaNonlinearOperator #
#####################################
"""
    struct GeneralizedAlphaNonlinearOperator <: DiscreteODEOperator end

Nonlinear discrete operator corresponding to the generalized-α scheme:
```math
residual(t_αf, u_n+αf, v_n+αm) = 0.
```
"""
struct GeneralizedAlphaNonlinearOperator <: DiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  vαm::AbstractVector
  t0::Real
  us0::NTuple{2,AbstractVector}
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  tαf::Real
end

function Algebra.allocate_residual(
  disop::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, dt, (u0, v0), vαm = disop.tαf, disop.dt, disop.us0, disop.vαm
  _uαf_to_vαm!(vαm, uαf, u0, v0, dt, disop.αm, disop.αf, disop.γ)
  allocate_residual(disop.odeop, t, (uαf, vαm), disop.odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  disop::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, dt, (u0, v0), vαm = disop.tαf, disop.dt, disop.us0, disop.vαm
  _uαf_to_vαm!(vαm, uαf, u0, v0, dt, disop.αm, disop.αf, disop.γ)
  residual!(r, disop.odeop, t, (uαf, vαm), disop.odeopcache)
end

function Algebra.allocate_jacobian(
  disop::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, dt, (u0, v0), vαm = disop.tαf, disop.dt, disop.us0, disop.vαm
  _uαf_to_vαm!(vαm, uαf, u0, v0, dt, disop.αm, disop.αf, disop.γ)
  allocate_jacobian(disop.odeop, t, (uαf, vαm), disop.odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  disop::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, dt, (u0, v0), vαm = disop.tαf, disop.dt, disop.us0, disop.vαm
  wu, vαm = _uαf_to_vαm!(vαm, uαf, u0, v0, dt, disop.αm, disop.αf, disop.γ)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, disop.odeop, t, (uαf, vαm), (1, wu), disop.odeopcache)
end

function Algebra.solve!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlphaNonlinearOperator,
  disslvrcache
)
  odeop, odeopcache = disop.odeop, disop.odeopcache
  t0, dt, (u0, v0) = disop.t0, disop.dt, disop.us0

  update_odeopcache!(odeopcache, odeop, t0)

  uF, vF = usF
  disslvrcache = solve!(uF, disslvr, disop, disslvrcache)
  _uF_to_usF!(uF, vF, u0, v0, dt, disop.αf, disop.γ)

  usF = (uF, vF)
  (usF, disslvrcache)
end

#####################################
# GeneralizedAlphaNonlinearOperator #
#####################################
"""
    struct GeneralizedAlphaLinearOperator <: LinearDiscreteODEOperator end

Linear discrete operator corresponding to the generalized-α scheme:
```math
residual(t_αf, u_n+αf, v_n+αm) = mass(t_αf) v_n+αm + stiffness(t_αf) u_n+αf + res(t_αf) = 0.
```
"""
struct GeneralizedAlphaLinearOperator <: LinearDiscreteODEOperator
  odeop::ODEOperator
  odeopcache
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  us0::NTuple{2,AbstractVector}
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  tαf::Real
end

Algebra.get_matrix(disop::GeneralizedAlphaLinearOperator) = disop.J

Algebra.get_vector(disop::GeneralizedAlphaLinearOperator) = disop.r

function Algebra.solve!(
  usF::NTuple{2,AbstractVector},
  disslvr::NonlinearSolver, disop::GeneralizedAlphaLinearOperator,
  disslvrcache
)
  odeop, odeopcache, J, r = disop.odeop, disop.odeopcache, disop.J, disop.r
  t0, dt, (u0, v0) = disop.t0, disop.dt, disop.us0
  αf, γ, tαf = disop.αf, disop.γ, disop.tαf

  update_odeopcache!(odeopcache, odeop, t0)

  wv, wu = 1 - disop.αm / disop.γ, disop.αm / (disop.γ * disop.αf * disop.dt)
  uF, vF = usF
  fill!(uF, zero(eltype(uF)))
  fill!(vF, zero(eltype(vF)))
  axpy!(wv, v0, vF)
  axpy!(-wu, u0, vF)

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, odeop, tαf, (uF, vF), (1, wu), odeopcache)
  residual!(r, odeop, tαf, (uF, vF), odeopcache)
  rmul!(r, -1)

  disslvrcache = solve!(uF, disslvr, disop, disslvrcache)
  _uF_to_usF!(uF, vF, u0, v0, dt, αf, γ)

  usF = (uF, vF)
  (usF, disslvrcache)
end

#########
# Utils #
#########
function _uF_to_usF!(uF, vF, u0, v0, dt, αf, γ)
  # @. uF = uF / αf + (1 - 1 / αf) * u0
  axpby!(1 - 1 / αf, u0, 1 / αf, uF)

  # @. vF = 1 / (γ * dt) * (uF - u0) + (1 - 1 / γ) * v0
  copy!(vF, uF)
  axpy!(-1, u0, vF)
  axpby!(1 - 1 / γ, v0, 1 / (γ * dt), vF)
  (uF, vF)
end

function _uαf_to_vαm!(vαm, uαf, u0, v0, dt, αm, αf, γ)
  # wv, wu = 1 - αm / γ, αm / (γ * αf * dt)
  # @. vαm = wv * v0 + wu * (uαf - u0)
  wv, wu = 1 - αm / γ, αm / (γ * αf * dt)
  copy!(vαm, uαf)
  axpy!(-1, u0, vαm)
  axpby!(wv, v0, wu, vαm)
  wu, vαm
end
