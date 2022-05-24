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

  (v0, ode_cache) = generalizedAlpha_cache
  ode_cache = update_cache!(ode_cache,op,tαf)
  nlop = GeneralizedAlphaNonlinearOperator(op,tαf,dt,αm,αf,γ,(u0,v0),generalizedAlpha_cache)
  nl_cache = solve!(u1,solver.nls,nlop,nl_cache)

  u1 = u1/αf + (1-1/αf)*u0
  v1 = 1/(γ*dt) * (u1-u0) + (1-1/γ)*v0

  cache = (generalizedAlpha_cache, nl_cache)
  x1 = (u1,v1)
  t1 = t0+dt

  return (x1,t1,cache)

end

"""
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
  vαm = (1 - op.αm/op.γ ) * v0 + op.αm/(op.γ*op.αf*op.dt) * (uαf - u0)
  residual!(b,op.odeop,op.tαf,(uαf,vαm),cache)
end

function jacobian!(A::AbstractMatrix,op::GeneralizedAlphaNonlinearOperator,x::AbstractVector)
  uαf = x
  u0, v0 = op.x0
  vαm, cache = op.ode_cache
  vαm = (1 - op.αm/op.γ ) * v0 + op.αm/(op.γ*op.αf*op.dt) * (uαf - u0)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.tαf,(uαf,vαm),(1.0,op.αm/(op.αf*op.γ*op.dt)),cache)
end

function allocate_residual(op::GeneralizedAlphaNonlinearOperator,x::AbstractVector)
  vαm, cache = op.ode_cache
  allocate_residual(op.odeop,x,cache)
end

function allocate_jacobian(op::GeneralizedAlphaNonlinearOperator,x::AbstractVector)
  vαm, cache = op.ode_cache
  allocate_jacobian(op.odeop,x,cache)
end

function zero_initial_guess(op::GeneralizedAlphaNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end
