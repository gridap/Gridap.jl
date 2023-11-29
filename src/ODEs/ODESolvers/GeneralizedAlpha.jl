"""
    struct GeneralizedAlpha <: ODESolver

Generalized-α ODE solver.
"""
struct GeneralizedAlpha <: ODESolver
  sol::NonlinearSolver
  dt::Float64
  ρ∞::Float64

  function GeneralizedAlpha(sol, dt, ρ∞)
    ρ∞01 = clamp(ρ∞, 0, 1)
    if ρ∞01 != ρ∞
      msg = """
      The parameter ρ∞ of the generalized-α must lie between zero and one.
      Setting ρ∞ to $(ρ∞01).
      """
      @warn msg
    end
    new(sol, dt, ρ∞01)
  end
end

# ODESolver interface
function get_dt(solver::GeneralizedAlpha)
  solver.dt
end

function allocate_dop_cache(
  solver::GeneralizedAlpha,
  ode_op::ODEOperator, ode_cache,
  t::Real, us::NTuple{2,AbstractVector}
)
  similar(us[2])
end

function allocate_dop_cache(
  solver::GeneralizedAlpha,
  ode_op::ODEOperator{LinearODE}, ode_cache,
  t::Real, us::NTuple{2,AbstractVector}
)
  J = allocate_jacobian(ode_op, t, us, ode_cache)
  r = allocate_residual(ode_op, t, us, ode_cache)
  (J, r)
end

function allocate_sol_cache(solver::GeneralizedAlpha)
  nothing
end

function DiscreteODEOperator(
  solver::GeneralizedAlpha, ode_op::ODEOperator, ode_cache,
  dop_cache, t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tαf::Real
)
  u0, v0 = us0
  v = dop_cache
  GeneralizedAlphaNonlinearOperator(
    ode_op, ode_cache, v,
    t0, u0, v0, dt, αm, αf, γ, tαf
  )
end

function DiscreteODEOperator(
  solver::GeneralizedAlpha, ode_op::ODEOperator{LinearODE}, ode_cache,
  dop_cache, t0::Real, us0::NTuple{2,AbstractVector}, dt::Real,
  αm::Real, αf::Real, γ::Real, tαf::Real
)
  u0, v0 = us0
  J, r = dop_cache
  GeneralizedAlphaLinearOperator(
    ode_op, ode_cache, J, r,
    t0, u0, v0, dt, αm, αf, γ, tαf
  )
end

###############
# solve_step! #
###############
function solve_step!(
  usF::NTuple{2,AbstractVector},
  solver::GeneralizedAlpha, ode_op::ODEOperator,
  us0::NTuple{2,AbstractVector}, t0::Real,
  cache
)
  # Unpack solver
  dt = get_dt(solver)
  ρ∞ = solver.ρ∞

  αf = 1 / (1 + ρ∞)
  αm = (3 - ρ∞) / (1 + ρ∞) / 2
  γ = 1 / 2 + αm - αf
  tαf = t0 + αf * dt

  u0, v0 = us0

  # Allocate or unpack cache
  if isnothing(cache)
    ode_cache = allocate_cache(ode_op, t0, (u0, v0))
    dop_cache = allocate_dop_cache(solver, ode_op, ode_cache, t0, us0)
    sol_cache = allocate_sol_cache(solver)
  else
    ode_cache, dop_cache, sol_cache = cache
  end

  # Create discrete ODE operator
  dop = DiscreteODEOperator(
    solver, ode_op, ode_cache, dop_cache,
    t0, us0, dt, αm, αf, γ, tαf
  )

  # Solve discrete ODE operator
  usF, sol_cache = solve_dop!(usF, dop, solver.sol, sol_cache)
  tF = t0 + dt

  # Update cache
  cache = (ode_cache, dop_cache, sol_cache)

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
  ode_op::ODEOperator
  ode_cache
  v::AbstractVector
  t0::Real
  u0::AbstractVector
  v0::AbstractVector
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  tαf::Real
end

function Algebra.allocate_residual(
  op::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, vαm = op.tαf, op.v
  _uαf_to_vαm!(vαm, uαf, op.u0, op.v0, op.dt, op.αm, op.αf, op.γ)
  allocate_residual(op.ode_op, t, (uαf, vαm), op.ode_cache)
end

function Algebra.residual!(
  r::AbstractVector,
  op::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, vαm = op.tαf, op.v
  _uαf_to_vαm!(vαm, uαf, op.u0, op.v0, op.dt, op.αm, op.αf, op.γ)
  residual!(r, op.ode_op, t, (uαf, vαm), op.ode_cache)
end

function Algebra.allocate_jacobian(
  op::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, vαm = op.tαf, op.v
  _uαf_to_vαm!(vαm, uαf, op.u0, op.v0, op.dt, op.αm, op.αf, op.γ)
  allocate_jacobian(op.ode_op, t, (uαf, vαm), op.ode_cache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  op::GeneralizedAlphaNonlinearOperator,
  u::AbstractVector
)
  uαf = u
  t, vαm = op.tαf, op.v
  wu, vαm = _uαf_to_vαm!(vαm, uαf, op.u0, op.v0, op.dt, op.αm, op.αf, op.γ)
  fillstored!(J, zero(eltype(J)))
  jacobians!(J, op.ode_op, t, (uαf, vαm), (1, wu), op.ode_cache)
end

function solve_dop!(
  usF::NTuple{2,AbstractVector},
  op::GeneralizedAlphaNonlinearOperator, sol::NonlinearSolver, cache
)
  ode_op, ode_cache = op.ode_op, op.ode_cache
  t0, u0, v0, dt = op.t0, op.u0, op.v0, op.dt
  αf, γ = op.αf, op.γ

  update_cache!(ode_cache, ode_op, t0)

  uF, vF = usF
  cache = solve!(uF, sol, op, cache)
  _uF_to_usF!(uF, vF, u0, v0, dt, αf, γ)

  usF = (uF, vF)
  (usF, cache)
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
  ode_op::ODEOperator
  ode_cache
  J::AbstractMatrix
  r::AbstractVector
  t0::Real
  u0::AbstractVector
  v0::AbstractVector
  dt::Real
  αm::Real
  αf::Real
  γ::Real
  tαf::Real
end

Algebra.get_matrix(op::GeneralizedAlphaLinearOperator) = op.J

Algebra.get_vector(op::GeneralizedAlphaLinearOperator) = op.r

function solve_dop!(
  usF::NTuple{2,AbstractVector},
  op::GeneralizedAlphaLinearOperator, sol::NonlinearSolver, cache
)
  ode_op, ode_cache, J, r = op.ode_op, op.ode_cache, op.J, op.r
  t0, u0, v0, dt = op.t0, op.u0, op.v0, op.dt
  αf, γ, tαf = op.αf, op.γ, op.tαf

  update_cache!(ode_cache, ode_op, t0)

  wv, wu = 1 - op.αm / op.γ, op.αm / (op.γ * op.αf * op.dt)
  uF, vF = usF
  fill!(uF, zero(eltype(uF)))
  fill!(vF, zero(eltype(vF)))
  axpy!(wv, v0, vF)
  axpy!(-wu, u0, vF)

  fillstored!(J, zero(eltype(J)))
  jacobians!(J, ode_op, tαf, (uF, vF), (1, wu), ode_cache)
  residual!(r, ode_op, tαf, (uF, vF), ode_cache)
  rmul!(r, -1)

  cache = solve!(uF, sol, op, cache)
  _uF_to_usF!(uF, vF, u0, v0, dt, αf, γ)

  usF = (uF, vF)
  (usF, cache)
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
  # wv, wu = 1 - op.αm / op.γ, op.αm / (op.γ * op.αf * op.dt)
  # @. vαm = wv * v0 + wu * (uαf - u0)
  wv, wu = 1 - αm / γ, αm / (γ * αf * dt)
  copy!(vαm, uαf)
  axpy!(-1, u0, vαm)
  axpby!(wv, v0, wu, vαm)
  wu, vαm
end
