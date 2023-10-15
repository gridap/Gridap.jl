
"""
Runge-Kutta ODE solver
"""
struct RungeKutta <: ODESolver
  nls_stage::NonlinearSolver
  nls_update::NonlinearSolver
  dt::Float64
  bt::ButcherTableau
  function RungeKutta(nls_stage::NonlinearSolver,nls_update::NonlinearSolver,dt,type::Symbol)
    bt = ButcherTableau(type)
    new(nls_stage,nls_update,dt,bt)
  end
end

function solve_step!(uf::AbstractVector,
  solver::RungeKutta,
  op::ODEOperator,
  u0::AbstractVector,
  t0::Real,
  cache)

  # Unpack variables
  dt = solver.dt
  s = solver.bt.s
  a = solver.bt.a
  b = solver.bt.b
  c = solver.bt.c
  d = solver.bt.d

  # Create cache if not there
  if cache === nothing
    ode_cache = allocate_cache(op)
    vi = similar(u0)
    fi = [similar(u0)]
    nls_stage_cache = nothing
    nls_update_cache = nothing
  else
    ode_cache, vi, fi, nls_stage_cache, nls_update_cache = cache
  end

  # Create RKNL stage operator
  nlop_stage = RungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,fi,0,a)

  # Compute intermediate stages
  for i in 1:s

    # allocate space to store the RHS at i
    if (length(fi) < i)
      push!(fi,similar(u0))
    end

    # Update time
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop_stage,ti,fi,i)

    if(a[i,i]==0)
      # Skip stage solve if a_ii=0 => u_i=u_0, f_i = f_0
      @. uf = u0
    else
      # solve at stage i
      nls_stage_cache = solve!(uf,solver.nls_stage,nlop_stage,nls_stage_cache)
    end

    # Update RHS at stage i using solution at u_i
    rhs!(nlop_stage, uf)

  end

  # Update final time
  tf = t0+dt

  # Skip final update if not necessary
  if !(c[s]==1.0 && a[s,:] == b)

    # Create RKNL final update operator
    ode_cache = update_cache!(ode_cache,op,tf)
    nlop_update = RungeKuttaUpdateNonlinearOperator(op,tf,dt,u0,ode_cache,vi,fi,s,b)

    # solve at final update
    nls_update_cache = solve!(uf,solver.nls_update,nlop_update,nls_update_cache)

  end

  # Update final cache
  cache = (ode_cache, vi, fi, nls_stage_cache, nls_update_cache)

  return (uf,tf,cache)

end

abstract type RungeKuttaNonlinearOperator <: NonlinearOperator end

"""
Nonlinear operator that represents the Runge-Kutta nonlinear operator at a
given time step and stage, i.e., A(t,u_i,(u_i-u_n)/dt)
"""
mutable struct RungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  fi::Vector{AbstractVector}
  i::Int
  a::Matrix
end

"""
Nonlinear operator that represents the Runge-Kutta nonlinear operator at the
final updated of a given time step, i.e., A(t,u_f,(u_f-u_n)/dt)
"""
mutable struct RungeKuttaUpdateNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  fi::Vector{AbstractVector}
  s::Int
  b::Vector{Number}
end



"""
Compute the residual of the Runge-Kutta nonlinear operator at stage i.
```math
A(t,ui,∂ui/∂t) = ∂ui/∂t - a_ii * f(ui,ti) - ∑_{j<i} a_ij * f(uj,tj) = 0
```

Uses the vector b as auxiliar variable to store the residual of the left-hand side of
the i-th stage ODE operator, then adds the corresponding contribution from right-hand side
at all earlier stages.
```math
b = M(ui,ti)∂u/∂t
b - ∑_{j<=i} a_ij * f(uj,tj) = 0
```
"""
function residual!(b::AbstractVector,op::RungeKuttaStageNonlinearOperator,x::AbstractVector)
  rhs!(op,x)
  lhs!(b,op,x)
  for j in 1:op.i
    @. b = b - op.a[op.i,j] * op.fi[j]
  end
  b
end

"""
Compute the residual of the final update of the Runge-Kutta nonlinear operator.
```math
A(t,uf,∂uf/∂t) = ∂uf/∂t - ∑_{i<=s} b_i * f(ui,ti) = 0
```

Uses the vector b as auxiliar variable to store the residual of the update ODE
operator (e.g. identity or mass matrix), then adds the corresponding contribution from earlier stages.
```math
b = [∂u/∂t]
b - ∑_{i<=s} bi * f(ui,ti) = 0
```
"""
function residual!(b::AbstractVector,op::RungeKuttaUpdateNonlinearOperator,x::AbstractVector)
  lhs!(b,op,x)
  for i in 1:op.s
    @. b = b - op.b[i] * op.fi[i]
  end
  b
end

function jacobian!(A::AbstractMatrix,op::RungeKuttaStageNonlinearOperator,x::AbstractVector)
  # @assert (abs(op.a[op.i,op.i]) > 0.0)
  ui = x
  vi = op.vi
  @. vi = (x-op.u0)/(op.dt)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.ti,(ui,vi),(op.a[op.i,op.i],1.0/op.dt),op.ode_cache)
end

function jacobian!(A::AbstractMatrix,op::RungeKuttaUpdateNonlinearOperator,x::AbstractVector)
  uf = x
  vf = op.vi
  @. vf = (x-op.u0)/(op.dt)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,op.odeop,op.ti,(uf,vf),2,1.0/(op.dt),op.ode_cache)
end

function jacobian!(A::AbstractMatrix,op::RungeKuttaStageNonlinearOperator,x::AbstractVector,
  i::Integer,γᵢ::Real)
  ui = x
  vi = op.vi
  @. vi = (x-op.u0)/(op.dt)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,op.odeop,op.ti,(ui,vi),i,γᵢ,op.ode_cache)
end

function allocate_residual(op::RungeKuttaNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function allocate_jacobian(op::RungeKuttaNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end

function zero_initial_guess(op::RungeKuttaNonlinearOperator)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end

function rhs!(op::RungeKuttaNonlinearOperator, x::AbstractVector)
  u = x
  v = op.vi
  @. v = (x-op.u0)/(op.dt)
  f = op.fi
  rhs!(f[op.i],op.odeop,op.ti,(u,v),op.ode_cache)
end

function lhs!(b::AbstractVector, op::RungeKuttaNonlinearOperator, x::AbstractVector)
  u = x
  v = op.vi
  @. v = (x-op.u0)/(op.dt)
  lhs!(b,op.odeop,op.ti,(u,v),op.ode_cache)
end

function update!(op::RungeKuttaNonlinearOperator,ti::Float64,fi::AbstractVector,i::Int)
  op.ti = ti
  op.fi = fi
  op.i = i
end

function _mass_matrix!(A,odeop,t,ode_cache,u)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,odeop,t,(u,u),2,1.0,ode_cache)
end
