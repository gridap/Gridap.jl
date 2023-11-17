
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
    ui = similar(u0)
    ki = Vector{typeof(u0)}()
    sizehint!(ki,s)
    [push!(ki,similar(u0)) for i in 1:s]
    rhs = similar(u0)
    nls_stage_cache = nothing
    nls_update_cache = nothing
  else
    ode_cache, ui, ki, rhs, nls_stage_cache, nls_update_cache = cache
  end

  # Initialize states to zero
  for i in 1:s
    @. ki[i] *= 0.0
  end

  # Create RKNL stage operator
  nlop_stage = RungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,ui,ki,rhs,0,a)

  # Compute intermediate stages
  for i in 1:s

    # Update time
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop_stage,ti,i)

    # solve at stage i
    nls_stage_cache = solve!(uf,solver.nls_stage,nlop_stage,nls_stage_cache)

    # Update stage unknown
    @. nlop_stage.ki[i] = uf

  end

  # Update final time
  tf = t0+dt

  # Update final solution
  @. uf = u0
  for i in 1:s
    @. uf = uf + dt * b[i] * nlop_stage.ki[i]
  end

  # Update final cache
  cache = (ode_cache, ui, ki, rhs, nls_stage_cache, nls_update_cache)

  return (uf,tf,cache)

end

abstract type RungeKuttaNonlinearOperator <: NonlinearOperator end

"""
Nonlinear operator that represents the Runge-Kutta nonlinear operator at a
given time step and stage, i.e., A(tᵢ,uᵢ,kᵢ)
"""
mutable struct RungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  ui::AbstractVector
  ki::Vector{AbstractVector}
  rhs::AbstractVector
  i::Int
  a::Matrix
end

"""
Compute the residual of the Runge-Kutta nonlinear operator at stage i.
```math
A(t,ui,ki) = M(ti) ki - f(u₀ + ∑_{j<=i} Δt * a_ij * kj, tj) = 0
```

Uses the vector b as auxiliar variable to store the residual of the left-hand side of
the i-th stage ODE operator, then adds the corresponding contribution from right-hand side
at all earlier stages.
```math
b = M(ti) Ki
b - f(u₀ + ∑_{j<=i} Δt * a_ij * kj, tj) = 0
```
"""
function residual!(b::AbstractVector,op::RungeKuttaStageNonlinearOperator,x::AbstractVector)
  rhs!(op,x)
  lhs!(b,op,x)
  @. b = b - op.rhs
  b
end

function jacobian!(A::AbstractMatrix,op::RungeKuttaStageNonlinearOperator,x::AbstractVector)
  u = op.ui
  @. u = op.u0
  for j in 1:op.i-1
    @. u = u + op.dt * op.a[op.i,j] * op.ki[j]
  end
  @. u = u + op.dt * op.a[op.i,op.i] * x
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.ti,(u,x),(op.dt*op.a[op.i,op.i],1.0),op.ode_cache)
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

function rhs!(op::RungeKuttaStageNonlinearOperator, x::AbstractVector)
  u = op.ui
  @. u = op.u0
  for j in 1:op.i-1
    @. u = u + op.dt * op.a[op.i,j] * op.ki[j]
  end
  @. u = u + op.dt * op.a[op.i,op.i] * x
  rhs!(op.rhs,op.odeop,op.ti,(u,x),op.ode_cache)
end

function lhs!(b::AbstractVector, op::RungeKuttaNonlinearOperator, x::AbstractVector)
  u = op.ui
  @. u *= 0
  lhs!(b,op.odeop,op.ti,(u,x),op.ode_cache)
end

function update!(op::RungeKuttaNonlinearOperator,ti::Float64,i::Int)
  op.ti = ti
  op.i = i
end

# Redefining solve! function to enforce computation of the jacobian within
# each stage of the Runge-Kutta method when the solver is "LinearSolver".
function solve!(x::AbstractVector,
  ls::LinearSolver,
  op::RungeKuttaNonlinearOperator,
  cache::Nothing)
  fill!(x,zero(eltype(x)))
  b = residual(op, x)
  A = jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)
  solve!(x,ns,b)
  LinearSolverCache(A,b,ns)
end

function solve!(x::AbstractVector,
  ls::LinearSolver,
  op::RungeKuttaNonlinearOperator,
  cache)
  fill!(x,zero(eltype(x)))
  b = cache.b
  A = cache.A
  ns = cache.ns
  residual!(b, op, x)
  jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)
  solve!(x,ns,b)
  cache
end
