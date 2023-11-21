"""
Generalized-α ODE solver
"""
struct GeneralizedAlpha <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  ρ∞::Float64
end

function solve_step!(
  x1::NTuple{2,AbstractVector},
  solver::GeneralizedAlpha,
  op::ODEOperator,
  x0::NTuple{2,AbstractVector},
  t0::Real,
  cache) # -> (uF,tF)

  dt = solver.dt
  ρ∞ = solver.ρ∞
  αf = 1.0/(1.0 + ρ∞)
  αm = 0.5 * (3-ρ∞) / (1+ρ∞)
  γ = 0.5 + αm - αf
  tαf = t0+αf*dt
  u0, v0 = x0
  u1, v1 = x1

  if cache === nothing
    generalizedAlpha_cache = allocate_cache(op,v0)
    nl_cache = nothing
  else
    generalizedAlpha_cache, nl_cache = cache
  end

  (v, ode_cache) = generalizedAlpha_cache
  ode_cache = update_cache!(ode_cache,op,tαf)
  nlop = GeneralizedAlphaNonlinearOperator(op,tαf,dt,αm,αf,γ,x0,generalizedAlpha_cache)
  nl_cache = solve!(u1,solver.nls,nlop,nl_cache)

  @. u1 = u1/αf + (1-1/αf)*u0
  @. v1 = 1/(γ*dt) * (u1-u0) + (1-1/γ)*v0

  cache = (generalizedAlpha_cache, nl_cache)
  x1 = (u1,v1)
  t1 = t0+dt

  return (x1,t1,cache)

end


function solve_step!(
  x1::NTuple{3,AbstractVector},
  solver::GeneralizedAlpha,
  op::ODEOperator,
  x0::NTuple{3,AbstractVector},
  t0::Real,
  cache) # -> (uF,tF)

  dt = solver.dt
  ρ∞ = solver.ρ∞
  αf = ρ∞/(ρ∞ + 1.0)
  αm = (2*ρ∞ - 1.0)/(ρ∞ + 1.0)
  γ = 0.5 - αm + αf
  β = 0.25*((1.0 - αm + αf)^2)
  tαf = t0 + (1.0-αf)*dt
  u0, v0, a0 = x0
  u1, v1, a1 = x1

  if cache === nothing
    generalizedAlphaDtt_cache = allocate_cache(op,v0,a0)
    nl_cache = nothing
  else
    generalizedAlphaDtt_cache, nl_cache = cache
  end

  (v, a, ode_cache) = generalizedAlphaDtt_cache
  ode_cache = update_cache!(ode_cache,op,tαf)
  nlop = GeneralizedAlphaDttNonlinearOperator(op,tαf,dt,αm,αf,γ,β,x0,generalizedAlphaDtt_cache)
  nl_cache = solve!(u1,solver.nls,nlop,nl_cache)


  @. u1 = 1.0 / (1.0 - αf) * u1 -
    αf / (1.0 - αf) * u0

  @. v1 = γ / (β*dt) * (u1 - u0) -
    (γ - β) / β * v0 -
    (γ - 2.0 * β) / (2.0 * β) * dt * a0

  @. a1 = 1.0 / (β * dt * dt) * (u1 - u0) -
    1.0 / (β * dt) * v0 -
    (1.0 - 2.0 * β) / (2.0 * β) * a0

  cache = (generalizedAlphaDtt_cache, nl_cache)
  x1 = (u1,v1,a1)
  t1 = t0+dt

  return (x1,t1,cache)

end


"""
Generalized-α 1st order ODE solver
Nonlinear operator that represents the Generalized-α method nonlinear operator at a
given time step, i.e., A(t_αf,u_n+αf,v_n+αm)
"""
struct GeneralizedAlphaNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  tαf::Float64
  dt::Float64
  αm::Float64
  αf::Float64
  γ::Float64
  x0::NTuple{2,AbstractVector}
  ode_cache
end

function residual!(b::AbstractVector,op::GeneralizedAlphaNonlinearOperator,x::AbstractVector)
  uαf = x
  u0, v0 = op.x0
  vαm, cache = op.ode_cache
  @. vαm = (1 - op.αm/op.γ ) * v0 + op.αm/(op.γ*op.αf*op.dt) * (uαf - u0)
  residual!(b,op.odeop,op.tαf,(uαf,vαm),cache)
end

function jacobian!(A::AbstractMatrix,op::GeneralizedAlphaNonlinearOperator,x::AbstractVector)
  uαf = x
  u0, v0 = op.x0
  vαm, cache = op.ode_cache
  @. vαm = (1 - op.αm/op.γ ) * v0 + op.αm/(op.γ*op.αf*op.dt) * (uαf - u0)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.tαf,(uαf,vαm),(1.0,op.αm/(op.αf*op.γ*op.dt)),cache)
end

function allocate_residual(op::GeneralizedAlphaNonlinearOperator,x::AbstractVector)
  vαm, cache = op.ode_cache
  allocate_residual(op.odeop,op.tαf,x,cache)
end

function allocate_jacobian(op::GeneralizedAlphaNonlinearOperator,x::AbstractVector)
  vαm, cache = op.ode_cache
  allocate_jacobian(op.odeop,op.tαf,x,cache)
end

function zero_initial_guess(op::GeneralizedAlphaNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end


"""
Generalized-α 2nd order ODE solver
Nonlinear operator that represents the Generalized-α method nonlinear operator at a
given time step, i.e., A(t_(n+1-αf),u_(n+1-αf),v_(n+1-αf),a_(n+1-αm))
"""

struct GeneralizedAlphaDttNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  tαf::Float64
  dt::Float64
  αm::Float64
  αf::Float64
  γ::Float64
  β::Float64
  x0::NTuple{3,AbstractVector}
  ode_cache
end


function residual!(b::AbstractVector,op::GeneralizedAlphaDttNonlinearOperator,x::AbstractVector)
  uαf = x
  u0, v0, a0 = op.x0
  vαf, aαm, cache = op.ode_cache

  @. vαf = (op.γ) / (op.β * op.dt) * (uαf - u0) +
    (op.αf - 1.0) * (op.γ - op.β) / op.β * v0 +
    (op.αf - 1.0) * (op.γ - 2.0*op.β) / (2.0 * op.β) * op.dt * a0 +
    op.αf * v0

  @. aαm = (1.0 - op.αm) / (1.0 - op.αf) / (op.β * op.dt * op.dt) * (uαf - u0) +
    (op.αm - 1.0) / (op.β * op.dt) * v0 +
    (op.αm - 1.0) * (1.0 - 2.0*op.β) / (2.0 * op.β) * a0 +
    op.αm * a0

  residual!(b,op.odeop,op.tαf,(uαf,vαf,aαm),cache)
end


function jacobian!(A::AbstractMatrix,op::GeneralizedAlphaDttNonlinearOperator,x::AbstractVector)
  uαf = x
  u0, v0, a0 = op.x0
  vαf, aαm, cache = op.ode_cache

  @. vαf = (op.γ) / (op.β * op.dt) * (uαf - u0) +
    (op.αf - 1.0) * (op.γ - op.β) / op.β * v0 +
    (op.αf - 1.0) * (op.γ - 2.0*op.β) / (2.0 * op.β) * op.dt * a0 +
    op.αf * v0

  @. aαm = (1.0 - op.αm) / (1.0 - op.αf) / (op.β * op.dt * op.dt) * (uαf - u0) +
    (op.αm - 1.0) / (op.β * op.dt) * v0 +
    (op.αm - 1.0) * (1.0 - 2.0*op.β) / (2.0 * op.β) * a0 +
    op.αm * a0

  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.tαf,(uαf,vαf,aαm),
    (1.0, op.γ/(op.β * op.dt),
      (1.0 - op.αm) / (1.0 - op.αf) / (op.β * op.dt * op.dt)),
    cache)
end


function allocate_residual(op::GeneralizedAlphaDttNonlinearOperator,x::AbstractVector)
  vαf, aαm, cache = op.ode_cache
  allocate_residual(op.odeop,op.tαf,x,cache)
end


function allocate_jacobian(op::GeneralizedAlphaDttNonlinearOperator,x::AbstractVector)
  vαf, aαm, cache = op.ode_cache
  allocate_jacobian(op.odeop,op.tαf,x,cache)
end


function zero_initial_guess(op::GeneralizedAlphaDttNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end
