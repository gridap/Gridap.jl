"""
Explicit Runge-Kutta ODE solver
"""
struct DIRungeKutta <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  tableau::ButcherTableau
  function DIRungeKutta(nls::NonlinearSolver, dt, type::Symbol)
    bt = ButcherTableau(type)
    new(nls, dt, bt)
  end
end

"""
solve_step!(uf,odesol,op,u0,t0,cache)
"""
function solve_step!(uf::AbstractVector,
  solver::DIRungeKutta,
  op::ODEOperator,
  u0::AbstractVector,
  t0::Real,
  cache)

  # Unpack variables
  dt = solver.dt
  s = solver.tableau.s
  a = solver.tableau.a
  b = solver.tableau.b
  c = solver.tableau.c
  d = solver.tableau.d

  # Create cache if not there
  if cache === nothing
    ode_cache = allocate_cache(op)
    vi = similar(u0)
    ki = [similar(u0) for i in 1:s]
    nl_cache = nothing
  else
    ode_cache, vi, ki, M, nl_cache = cache
  end

  nlop = DIRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,ki,0,a)

  for i in 1:s

    # solve at stage i
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop,ti,ki[i],i)
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)

    update!(nlop,ti,uf,i)

  end

  # update final solution
  tf = t0 + dt

  @. uf = u0
  for i in 1:s
  @. uf = uf + dt*b[i]*nlop.ki[i]
  end

  cache = (ode_cache, vi, ki, M, nl_cache)

  return (uf,tf,cache)


end




mutable struct DIRungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  ki::AbstractVector
  i::Int
  a::Matrix
end


"""
ODE: A(t,u,∂u) = M ∂u/∂t + K(t,u) = 0 -> solve for u
RK:  A(t,u,ki) = M ki  + K(ti, dt a_ij ki)  + K(ti,u0 + dt ∑_{j<i} a_ij * kj) = 0 -> solve for ki
               = M ki    + K(ti,ui) = 0
For forward euler, i = 1     -> ui = u0
For other methods, i = 1,…,s -> ui = u0 + dt ∑_{j<i} a_ij * kj
"""
function residual!(b::AbstractVector,op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)

  ui = x
  vi = op.vi

  lhs!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. ui = op.u0
  for j = 1:op.i-1
   @. ui = ui  + op.dt * op.a[op.i,j] * op.ki[j]
  end

  @. ui = ui + op.dt * op.a[op.i,i] * x

  rhs = similar(op.u0)
  rhs!(rhs,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. b = b + rhs
  @. b = -1.0 * b
  b
end

function jacobian!(A::AbstractMatrix,op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)

  ui = x
  vi = op.vi
  @. ui = x # this value is irrelevant its jacobian contribution is zero
  @. vi = x
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.ti,(ui,vi),(op.dt*op.a[op.i,i],1.0),op.ode_cache)

end


function allocate_residual(op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function allocate_jacobian(op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end


function update!(op::DIRungeKuttaStageNonlinearOperator,ti::Float64,ki::AbstractVector,i::Int)
  op.ti = ti
  @. op.ki[i] = ki
  op.i = i
end
