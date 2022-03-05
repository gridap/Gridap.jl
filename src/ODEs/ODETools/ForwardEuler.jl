"""
f-method ODE solver
"""
struct ForwardEuler <: ODESolver
  nls::NonlinearSolver
  dt::Float64
end

function solve_step!(uf::AbstractVector,
                     solver::ForwardEuler,
                     op::ODEOperator,
                     u0::AbstractVector,
                     t0::Real,
                     cache) # -> (uF,tF)

  if cache === nothing
    ode_cache = allocate_cache(op)
    vf = similar(u0)
    nl_cache = nothing
  else
    ode_cache, vf, nl_cache = cache
  end

  dt = solver.dt
  tf = t0+dt
  # The space should have the boundary conditions at tf
  ode_cache = update_cache!(ode_cache,op,t0)

  nlop = ForwardEulerNonlinearOperator(op,t0,dt,u0,ode_cache,vf)

  nl_cache = solve!(uf,solver.nls,nlop,nl_cache)

  cache = (ode_cache, vf, nl_cache)

  return (uf,tf,cache)

end

"""
Nonlinear operator that represents the Forward Euler nonlinear operator at a
given time step, i.e., A(t,u_n,(u_n+1-u_n)/dt)
"""
struct ForwardEulerNonlinearOperator <: NonlinearOperator
  odeop::ODEOperator
  tf::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vf::AbstractVector
end

function residual!(b::AbstractVector,op::ForwardEulerNonlinearOperator,x::AbstractVector)
  vf = op.vf
  vf = (x-op.u0)/op.dt
  residual!(b,op.odeop,op.tf,(op.u0,vf),op.ode_cache)
end

function jacobian!(A::AbstractMatrix,op::ForwardEulerNonlinearOperator,x::AbstractVector)
  vf = op.vf
  vf = (x-op.u0)/op.dt
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.tf,(op.u0,vf),(0,1/op.dt),op.ode_cache)
end

function allocate_residual(op::ForwardEulerNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,x,op.ode_cache)
end

function allocate_jacobian(op::ForwardEulerNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,x,op.ode_cache)
end

function zero_initial_guess(op::ForwardEulerNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end
