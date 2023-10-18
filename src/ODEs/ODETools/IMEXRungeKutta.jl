"""
Implicit-Explicit Runge-Kutta ODE solver.

This struct defines an ODE solver for the system of ODEs

    M(u,t)du/dt = f(u,t) + g(u,t)

where `f` is a nonlinear function of `u` and `t` that will treated implicitly and
  `g` is a nonlinear function of `u` and `t` that will be treated explicitly.
  The ODE is solved using an implicit-explicit Runge-Kutta method.
"""
struct IMEXRungeKutta <: ODESolver
  nls_stage::NonlinearSolver
  nls_update::NonlinearSolver
  dt::Float64
  tableau::IMEXButcherTableau
  function IMEXRungeKutta(nls_stage::NonlinearSolver, nls_update::NonlinearSolver, dt, type::Symbol)
    bt = IMEXButcherTableau(type)
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
  solver::IMEXRungeKutta,
  op::ODEOperator,
  u0::AbstractVector,
  t0::Real,
  cache)

  # Unpack variables
  dt = solver.dt
  s = solver.tableau.s
  aᵢ = solver.tableau.aᵢ
  bᵢ = solver.tableau.bᵢ
  aₑ = solver.tableau.aₑ
  bₑ = solver.tableau.bₑ
  c = solver.tableau.c
  d = solver.tableau.d

  # Create cache if not there
  if cache === nothing
    ode_cache = allocate_cache(op)
    vi = similar(u0)
    ui = Vector{typeof(u0)}(undef,0)
    im_rhs = similar(u0)
    ex_rhs = similar(u0)
    for i in 1:s
      push!(ui,similar(u0))
    end
    nls_stage_cache = nothing
    nls_update_cache = nothing
  else
    ode_cache, vi, ui, im_rhs, ex_rhs, nls_stage_cache, nls_update_cache = cache
  end

  # Create RKNL stage operator
  nlop_stage = IMEXRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,ui,im_rhs,ex_rhs,0,aᵢ,aₑ)

  # Compute intermediate stages
  for i in 1:s

    # Update time
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop_stage,ti,ui,i)

    if(aᵢ[i,i]==0)
      # Skip stage solve if a_ii=0 => u_i=u_0, f_i = f_0, gi = g_0
      @. uf = u0
    else
      # solve at stage i
      nls_stage_cache = solve!(uf,solver.nls_stage,nlop_stage,nls_stage_cache)
    end

    # Update stage unknown
    @. nlop_stage.ui[i] = uf

  end

  # Update final time
  tf = t0+dt

  # Skip final update if not necessary
  if !(c[s]==1.0 && aᵢ[s,:] == bᵢ && aₑ[s,:] == bₑ)

    # Create RKNL final update operator
    ode_cache = update_cache!(ode_cache,op,tf)
    nlop_update = IMEXRungeKuttaUpdateNonlinearOperator(op,tf,dt,u0,ode_cache,vi,ui,im_rhs,ex_rhs,s,bᵢ,bₑ)

    # solve at final update
    nls_update_cache = solve!(uf,solver.nls_update,nlop_update,nls_update_cache)

  end

  # Update final cache
  cache = (ode_cache, vi, ui, im_rhs, ex_rhs, nls_stage_cache, nls_update_cache)

  return (uf, tf, cache)

end

"""
IMEXRungeKuttaStageNonlinearOperator <: NonlinearOperator

Nonlinear operator for the implicit-explicit Runge-Kutta stage.
  At a given stage `i` it represents the nonlinear operator A(t,u_i,(u_i-u_n)/dt) such that
```math
A(t,u_i,(u_i-u_n)/dt) = M(u_i,t)(u_i-u_n)/Δt - f(∑aᵢ[i,j] * u_j,t_i) - g(∑aₑ[i,j] * u_j,t_i) = 0
```
"""
mutable struct IMEXRungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  ui::Vector{AbstractVector}
  im_rhs::AbstractVector
  ex_rhs::AbstractVector
  i::Int
  aᵢ::Matrix{Float64}
  aₑ::Matrix{Float64}
end

"""
IMEXRungeKuttaUpdateNonlinearOperator <: NonlinearOperator

Nonlinear operator for the implicit-explicit Runge-Kutta final update.
  At the final update it represents the nonlinear operator A(t,u_t,(u_t-u_n)/dt)  such that
  ```math
  A(t,u_f,(u_f-u_n)/dt) = M(u_f,t)(u_f-u_n)/Δt - ∑aᵢ[i,j] * f(u_j,t_j) - ∑aₑ[i,j] * g(u_j,t_j) = 0
  ```
"""
mutable struct IMEXRungeKuttaUpdateNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  ui::Vector{AbstractVector}
  im_rhs::AbstractVector
  ex_rhs::AbstractVector
  s::Int
  bᵢ::Vector{Float64}
  bₑ::Vector{Float64}
end

"""
residual!(b,op::IMEXRungeKuttaStageNonlinearOperator,x)

Compute the residual of the IMEXR Runge-Kutta nonlinear operator `op` at `x` and
store it in `b` for a given stage `i`.
```math
b = A(t,x,(x-x₀)/dt) = M(ui,ti)∂ui/∂t - f(∑aᵢ[i,j] * xj,ti) - g(∑_{j<i} aₑ_ij * xj,ti)
```

Uses the vector b as auxiliar variable to store the residual of the left-hand side of
the i-th stage ODE operator, then adds the corresponding contribution from right-hand side
at all earlier stages.
```math
b = M(ui,ti)∂u/∂t
b - f(∑_{j<=i} aᵢ_ij * uj,ti) - g(∑_{j<i} aₑ_ij * uj,ti) = 0
```
"""
function residual!(b::AbstractVector,
  op::IMEXRungeKuttaStageNonlinearOperator,
  x::AbstractVector)
  rhs!(op,x)
  explicit_rhs!(op,x)
  lhs!(b,op,x)
  @. b = b - op.im_rhs - op.ex_rhs
  b
end

"""
residual!(b,op::IMEXRungeKuttaUpdateNonlinearOperator,x)

Computes the residual of the IMEX Runge-Kutta nonlinear operator `op` at `x` and
for the final update.
```math
b = A(t,x,(x-x₀)/dt) = M(xf,tf)∂uf/∂t - f(∑bᵢ[i] * xi,tf) - g(∑bₑ[i] * xi,tf)
```

Uses the vector b as auxiliar variable to store the residual of the left-hand side of
the final update ODE operator, then adds the corresponding contribution from right-hand sides
at all stages.
```math
b = M(uf,tf)∂u/∂t
b - f(∑_{i<=s} bᵢ[i] * ui,tf) - g(∑_{i<=s} bₑ[i] * ui,tf) = 0
```
"""
function residual!(b::AbstractVector,
  op::IMEXRungeKuttaUpdateNonlinearOperator,
  x::AbstractVector)
  rhs!(op,x)
  explicit_rhs!(op,x)
  lhs!(b,op,x)
  @. b = b - op.im_rhs - op.ex_rhs
  b
end

"""
jacobian!(J,op::IMEXRungeKuttaStageNonlinearOperator,x)

Computes the Jacobian of the IMEX Runge-Kutta nonlinear operator `op` at `x` and
stores it in `J` for a given stage `i`.
```math
J = ∂A(t,x,(x-x₀)/dt)/∂x = ∂(xi/∂t)/∂xi + aᵢ[i,i] * (- ∂f(xi,ti)/∂xi )
```

It uses the function `jacobians!` to compute the Jacobian of the right-hand side (stiffness matrix)
and the jacobian of the left-hand side (mass matrix) at the current stage added together with the
coefficients `aᵢ[i,i]` and `1/Δt`, respectively. It assumes that the jacobian of the right-hand side
has already the negative sign incorporated in its calculation.
"""
function jacobian!(J::AbstractMatrix,
  op::IMEXRungeKuttaStageNonlinearOperator,
  x::AbstractVector)
  ui = x
  vi = op.vi
  @. vi = (x-op.u0)/(op.dt)
  z = zero(eltype(J))
  fillstored!(J,z)
  jacobians!(J,op.odeop,op.ti,(ui,vi),(op.aᵢ[op.i,op.i],1.0/op.dt),op.ode_cache)
end

"""
jacobian!(J,op::IMEXRungeKuttaUpdateNonlinearOperator,x)

Computes the Jacobian of the IMEX Runge-Kutta nonlinear operator `op` at `x` and
stores it in `J` for the final update. This is typically the mass matrix
"""
function jacobian!(J::AbstractMatrix,
  op::IMEXRungeKuttaUpdateNonlinearOperator,
  x::AbstractVector)
  uf = x
  vf = op.vi
  @. vf = (x-op.u0)/(op.dt)
  z = zero(eltype(J))
  fillstored!(J,z)
  jacobian!(J,op.odeop,op.ti,(uf,vf),2,1.0/(op.dt),op.ode_cache)
end

function rhs!(op::IMEXRungeKuttaStageNonlinearOperator, x::AbstractVector)
  u = zeros(eltype(x),length(x))
  for j in 1:op.i
    @. u = u + op.aᵢ[op.i,j] * op.ui[j]
  end
  v = op.vi
  @. v = (x-op.u0)/(op.dt)
  rhs!(op.im_rhs,op.odeop,op.ti,(u,v),op.ode_cache)
end

function rhs!(op::IMEXRungeKuttaUpdateNonlinearOperator, x::AbstractVector)
  u = zeros(eltype(x),length(x))
  for i in 1:op.s
    @. u = u + op.bᵢ[i] * op.ui[i]
  end
  v = op.vi
  @. v = (x-op.u0)/(op.dt)
  rhs!(op.im_rhs,op.odeop,op.ti,(u,v),op.ode_cache)
end

function explicit_rhs!(op::IMEXRungeKuttaStageNonlinearOperator, x::AbstractVector)
  u = zeros(eltype(x),length(x))
  for j in 1:op.i-1
    @. u = u + op.aₑ[op.i,j] * op.ui[j]
  end
  v = op.vi
  @. v = (x-op.u0)/(op.dt)
  explicit_rhs!(op.ex_rhs,op.odeop,op.ti,(u,v),op.ode_cache)
end

function explicit_rhs!(op::IMEXRungeKuttaUpdateNonlinearOperator, x::AbstractVector)
  u = zeros(eltype(x),length(x))
  for i in 1:op.s
    @. u = u + op.bₑ[i] * op.ui[i]
  end
  v = op.vi
  @. v = (x-op.u0)/(op.dt)
  explicit_rhs!(op.ex_rhs,op.odeop,op.ti,(u,v),op.ode_cache)
end
