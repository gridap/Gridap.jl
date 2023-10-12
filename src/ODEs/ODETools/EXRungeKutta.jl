"""
Explicit Runge-Kutta ODE solver.

This struct defines an ODE solver for the system of ODEs

    M(u,t)du/dt = g(u,t)

where `g` is a nonlinear function of `u` and `t` that will be treated explicitly.
  The ODE is solved using an explicit Runge-Kutta method.
"""
struct EXRungeKutta <: ODESolver
  nls_stage::NonlinearSolver
  nls_update::NonlinearSolver
  dt::Float64
  tableau::EXButcherTableau
  function EXRungeKutta(nls_stage::NonlinearSolver, nls_update::NonlinearSolver, dt, type::Symbol)
    bt = EXButcherTableau(type)
    new(nls_stage, nls_update, dt, bt)
  end
end

"""
solve_step!(uf,odesol,op,u0,t0,cache)

Solve one step of the ODE problem defined by `op` using the ODE solver `odesol`
  with initial solution `u0` at time `t0`. The solution is stored in `uf` and
  the final time in `tf`. The cache is used to store the solution of the
  nonlinear system of equations and auxiliar variables.
"""
function solve_step!(uf::AbstractVector,
  solver::EXRungeKutta,
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
    fi = Vector{typeof(u0)}(undef,0)
    for i in 1:s
      push!(fi,similar(u0))
    end
    nls_stage_cache = nothing
    nls_update_cache = nothing
  else
    ode_cache, vi, fi, nls_stage_cache, nls_update_cache = cache
  end

  # Create RKNL stage operator
  nlop_stage = EXRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,fi,0,a)

  # Compute intermediate stages
  for i in 1:s

    # Update time
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop_stage,ti,fi,i)

    if(i==0)
      # First stage is always forward euler: u_i = u_0 + dt*f_0
      @. uf = u0 + dt*fi
    else
      # solve at stage i
      nls_stage_cache = solve!(uf,solver.nls_stage,nlop_stage,nls_stage_cache)
    end

    # Update RHS at stage i using solution at u_i
    rhs!(nlop_stage, uf)
    # explicit_rhs!(nlop_stage, uf)

  end

  # Update final time
  tf = t0+dt

  # Skip final update if not necessary
  if !(c[s]==1.0 && a[s,:] == b)

    # Create RKNL final update operator
    ode_cache = update_cache!(ode_cache,op,tf)
    nlop_update = EXRungeKuttaUpdateNonlinearOperator(op,tf,dt,u0,ode_cache,vi,fi,s,b)

    # solve at final update
    nls_update_cache = solve!(uf,solver.nls_update,nlop_update,nls_update_cache)

  end

  # Update final cache
  cache = (ode_cache, vi, fi, nls_stage_cache, nls_update_cache)

  return (uf, tf, cache)

end

"""
EXRungeKuttaStageNonlinearOperator <: NonlinearOperator

Nonlinear operator for the explicit Runge-Kutta stage.
  At a given stage `i` it represents the nonlinear operator A(t,u_i,(u_i-u_n)/dt) such that
```math
A(t,u_i,(u_i-u_n)/dt) = M(u_i,t)(u_i-u_n)/Δt - ∑a[i,j] * g(u_j,t_j) = 0
```
"""
mutable struct EXRungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  fi::Vector{AbstractVector}
  i::Int
  a::Matrix{Float64}
end

"""
EXRungeKuttaUpdateNonlinearOperator <: NonlinearOperator

Nonlinear operator for the explicit Runge-Kutta final update.
  At the final update it represents the nonlinear operator A(t,u_t,(u_t-u_n)/dt)  such that
  ```math
  A(t,u_f,(u_f-u_n)/dt) = M(u_f,t)(u_f-u_n)/Δt - ∑aₑ[i,j] * g(u_j,t_j) = 0
  ```
"""
mutable struct EXRungeKuttaUpdateNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  fi::Vector{AbstractVector}
  s::Int
  b::Vector{Float64}
end

"""
residual!(b,op::IMEXRungeKuttaStageNonlinearOperator,x)

Compute the residual of the EX Runge-Kutta nonlinear operator `op` at `x` and
store it in `b` for a given stage `i`.
```math
b = A(t,x,(x-x₀)/dt) = ∂ui/∂t - ∑aᵢ[i,j] * f(xj,tj)
```

Uses the vector b as auxiliar variable to store the residual of the left-hand side of
the i-th stage ODE operator, then adds the corresponding contribution from right-hand side
at all earlier stages.
```math
b = M(ui,ti)∂u/∂t
b -  ∑_{j<i} aₑ_ij * g(uj,tj) = 0
```
"""
function residual!(b::AbstractVector,
  op::EXRungeKuttaStageNonlinearOperator,
  x::AbstractVector)
  rhs!(op,x)
  lhs!(b,op,x)
  for j in 1:op.i
    @. b = b - op.a[op.i,j] * op.fi[j]
  end
  b
end

"""
residual!(b,op::EXRungeKuttaUpdateNonlinearOperator,x)

Computes the residual of the EX Runge-Kutta nonlinear operator `op` at `x` and
for the final update.
```math
b = A(t,x,(x-x₀)/dt) = ∂ui/∂t - ∑bₑ[i] * g(xj,tj)
```

Uses the vector b as auxiliar variable to store the residual of the left-hand side of
the final update ODE operator, then adds the corresponding contribution from right-hand sides
at all stages.
```math
b = M(ui,ti)∂u/∂t
b - ∑_{i<s} bₑ[i] * g(ui,ti) = 0
```
"""
function residual!(b::AbstractVector,
  op::EXRungeKuttaUpdateNonlinearOperator,
  x::AbstractVector)
  lhs!(b,op,x)
  for i in 1:op.s
    @. b = b - op.b[i] * op.fi[i]
  end
  b
end

"""
jacobian!(J,op::EXRungeKuttaStageNonlinearOperator,x)

Computes the Jacobian of the EX Runge-Kutta nonlinear operator `op` at `x` and
stores it in `J` for a given stage `i`.
```math
J = ∂A(t,x,(x-x₀)/dt)/∂x = ∂(xi/∂t)/∂xi
```

It uses the function `jacobians!` to compute the Jacobian of the right-hand side (stiffness matrix)
and the jacobian of the left-hand side (mass matrix) at the current stage added together with the
coefficients `aᵢ[i,i]` and `1/Δt`, respectively. It assumes that the jacobian of the right-hand side
has already the negative sign incorporated in its calculation.
"""
function jacobian!(J::AbstractMatrix,
  op::EXRungeKuttaStageNonlinearOperator,
  x::AbstractVector)
  ui = x
  vi = op.vi
  @. vi = (x-op.u0)/(op.dt)
  z = zero(eltype(J))
  fillstored!(J,z)
  jacobians!(J,op.odeop,op.ti,(ui,vi),(op.aᵢ[op.i,op.i],1.0/op.dt),op.ode_cache)
end

"""
jacobian!(J,op::EXRungeKuttaUpdateNonlinearOperator,x)

Computes the Jacobian of the EX Runge-Kutta nonlinear operator `op` at `x` and
stores it in `J` for the final update. This is typically the mass matrix
"""
function jacobian!(J::AbstractMatrix,
  op::EXRungeKuttaUpdateNonlinearOperator,
  x::AbstractVector)
  uf = x
  vf = op.vi
  @. vf = (x-op.u0)/(op.dt)
  z = zero(eltype(J))
  fillstored!(J,z)
  jacobian!(J,op.odeop,op.ti,(uf,vf),2,1.0/(op.dt),op.ode_cache)
end

# function explicit_rhs!(op::RungeKuttaNonlinearOperator, x::AbstractVector)
#   u = x
#   v = op.vi
#   @. v = (x-op.u0)/(op.dt)
#   g = op.fi
#   explicit_rhs!(g[op.i],op.odeop,op.ti,(u,v),op.ode_cache)
# end

function update!(op::RungeKuttaNonlinearOperator,ti::Float64,fi::AbstractVector,i::Int)
  op.ti = ti
  op.fi = fi
  op.i = i
end
