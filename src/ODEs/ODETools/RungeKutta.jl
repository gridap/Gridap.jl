abstract type ButcherTableauType end

struct BE_1_0_1 <: ButcherTableauType end
struct CN_2_0_2 <: ButcherTableauType end
struct SDIRK_2_0_2 <: ButcherTableauType end
struct SDIRK_2_0_3 <: ButcherTableauType end
struct ESDIRK_3_1_2 <: ButcherTableauType end
struct TRBDF2_3_2_3 <: ButcherTableauType end

"""
Butcher tableau
"""
struct ButcherTableau{T <: ButcherTableauType}
  s::Int # stages
  p::Int # embedded order
  q::Int # order
  a::Matrix # A_ij
  b::Vector # b_j
  c::Vector # c_i
  d::Vector # d_j (embedded)
end

# Butcher Tableaus constructors
"""
Backward-Euler

number of stages: 1
embedded method: no
order: 1
"""
function ButcherTableau(::BE_1_0_1)
  s = 1
  p = 0
  q = 1
  a = reshape([1.0],1,1)
  b = [1.0]
  c = [1.0]
  d = [0.0]
  ButcherTableau{BE_1_0_1}(s,p,q,a,b,c,d)
end

"""
Crank-Nicolson (equivalent to trapezoidal rule)

number of stages: 2
embedded method: no
order: 2
"""
function ButcherTableau(type::CN_2_0_2)
  s = 2
  p = 0
  q = 2
  a = [0.0 0.0; 0.5 0.5]
  b = [0.5, 0.5]
  c = [0.0, 1.0]
  d = [0.0, 0.0]
  ButcherTableau{CN_2_0_2}(s,p,q,a,b,c,d)
end

"""
Qin and Zhang's SDIRK

number of stages: 2
embedded method: no
order: 2
"""
function ButcherTableau(type::SDIRK_2_0_2)
  s = 2
  p = 0
  q = 2
  a = [0.25 0.0; 0.5 0.25]
  b = [0.5, 0.5]
  c = [0.25, 0.75]
  d = [0.0, 0.0]
  ButcherTableau{SDIRK_2_0_2}(s,p,q,a,b,c,d)
end

"""
3rd order SDIRK

number of stages: 2
embedded method: no
order: 3
"""
function ButcherTableau(type::SDIRK_2_0_3)
  s = 2
  p = 0
  q = 3
  γ = (3-√(3))/6
  a = [γ 0.0; 1-2γ γ]
  b = [0.5, 0.5]
  c = [γ, 1-γ]
  d = [0.0, 0.0]
  ButcherTableau{SDIRK_2_0_3}(s,p,q,a,b,c,d)
end

function ButcherTableau(type::ESDIRK_3_1_2)
s = 3
p = 1
q = 2
γ = (2-√(2))/2
b₂ = (1 − 2γ)/(4γ)
b̂₂ = γ*(−2 + 7γ − 5(γ^2) + 4(γ^3)) / (2(2γ − 1))
b̂₃ = −2*(γ^2)*(1 − γ + γ^2) / (2γ − 1)
a = [0.0 0.0 0.0; γ γ 0.0; (1 − b₂ − γ) b₂ γ]
b = [(1 − b₂ − γ), b₂, γ]
c = [0.0, 2γ, 1.0]
d = [(1 − b̂₂ − b̂₃), b̂₂, b̂₃]
ButcherTableau{ESDIRK_3_1_2}(s,p,q,a,b,c,d)
end

function ButcherTableau(type::TRBDF2_3_2_3)
  s = 3
  p = 2
  q = 3
  aux = 2.0-√2.0
  a = [0.0 0.0 0.0; aux/2 aux/2 0.0; √2/4 √2/4 aux/2]
  b = [√2/4, √2/4, aux/2]
  c = [0.0, aux, 1.0]
  d = [(1.0-(√2/4))/3, ((3*√2)/4+1.0)/3, aux/6]
  ButcherTableau{TRBDF2_3_2_3}(s,p,q,a,b,c,d)
end

function ButcherTableau(type::Symbol)
  eval(:(ButcherTableau($type())))
end

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
    # M = allocate_jacobian(op,t0,u0,ode_cache)
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
  fi::AbstractVector
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

Uses the vector b as auxiliar variable to store the residual of the i-th stage
ODE operator, then adds the corresponding contribution from earlier stages.
```math
b = [1/a_ii * ∂u/∂t - f(ui,ti)]
Res_ij = - a_ij/a_ii * f(uj,ti)
b + ∑_{j<i} Res_ij = 0
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
b - ∑_{i<s} bi * f(ui,ti) = 0
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
