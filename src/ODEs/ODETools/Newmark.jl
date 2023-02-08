"""
Newmark-beta 2nd order ODE solver
"""
struct Newmark <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  γ::Float64
  β::Float64
end

function solve_step!(
  x1::NTuple{3,AbstractVector},
  solver::Newmark,
  op::ODEOperator,
  x0::NTuple{3,AbstractVector},
  t0::Real,
  cache) # -> (uF,tF)

  dt = solver.dt
  γ = solver.γ
  β = solver.β
  t1 = t0+dt
  u0, v0, a0 = x0
  u1, v1, a1 = x1

  if cache === nothing
    newmark_cache = allocate_cache(op,v0,a0)
    nl_cache = nothing
  else
    newmark_cache, nl_cache = cache
  end

  (v,a, ode_cache) = newmark_cache
  ode_cache = update_cache!(ode_cache,op,t1)
  nlop = NewmarkNonlinearOperator(op,t1,dt,γ,β,(u0,v0,a0),newmark_cache)
  nl_cache = solve!(u1,solver.nls,nlop,nl_cache)

  v1 = γ/(β*dt)*(u1-u0) + (1-γ/β)*v0 + dt*(1-γ/(2*β))*a0
  a1 = 1.0/(β*dt^2)*(u1-u0) - 1.0/(β*dt)*v0 - (1-2*β)/(2*β)*a0

  cache = (newmark_cache, nl_cache)
  x1 = (u1,v1,a1)

  return (x1,t1,cache)

end

"""
Nonlinear operator that represents the Newmark nonlinear operator at a
given time step, i.e., A(t,u_n+1,v_n+1,a_n+1)
"""
struct NewmarkNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  t1::Float64
  dt::Float64
  γ::Float64
  β::Float64
  x0::NTuple{3,AbstractVector}
  ode_cache
end

function residual!(b::AbstractVector,op::NewmarkNonlinearOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  residual!(b,op.odeop,op.t1,(u1,v1,a1),cache)
end

function jacobian!(A::AbstractMatrix,op::NewmarkNonlinearOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.t1,(u1,v1,a1),(1.0,op.γ/(op.β*op.dt),1.0/(op.β*op.dt^2)),cache)
end

function allocate_residual(op::NewmarkNonlinearOperator,x::AbstractVector)
  v1, a1, cache = op.ode_cache
  allocate_residual(op.odeop,x,cache)
end

function allocate_jacobian(op::NewmarkNonlinearOperator,x::AbstractVector)
  v1, a1, cache = op.ode_cache
  allocate_jacobian(op.odeop,x,cache)
end

function zero_initial_guess(op::NewmarkNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end
