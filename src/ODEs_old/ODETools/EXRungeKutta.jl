"""
Explicit Runge-Kutta ODE solver
"""
struct EXRungeKutta <: ODESolver
  ls::LinearSolver
  dt::Float64
  tableau::EXButcherTableau
  function EXRungeKutta(ls::LinearSolver, dt, type::Symbol)
    bt = EXButcherTableau(type)
    new(ls, dt, bt)
  end
end

"""
solve_step!(uf,odesol,op,u0,t0,cache)
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
    ki = [similar(u0) for i in 1:s]
    M = allocate_jacobian(op,t0,uf,ode_cache)
    get_mass_matrix!(M,op,t0,uf,ode_cache)
    l_cache = nothing
  else
    ode_cache, vi, ki, M, l_cache = cache
  end

  lop = EXRungeKuttaStageOperator(op,t0,dt,u0,ode_cache,vi,ki,0,a,M)

  for i in 1:s

    # solve at stage i
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(lop,ti,ki[i],i)
    l_cache = solve!(uf,solver.ls,lop,l_cache)

    update!(lop,ti,uf,i)

  end

  # update final solution
  tf = t0 + dt

  @. uf = u0
  for i in 1:s
  @. uf = uf + dt*b[i]*lop.ki[i]
  end

  cache = (ode_cache, vi, ki, M, l_cache)

  return (uf,tf,cache)


end




mutable struct EXRungeKuttaStageOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
  ki::AbstractVector
  i::Int
  a::Matrix
  M::AbstractMatrix
end


"""
ODE:    A(t,u,∂u  = M ∂u/∂t + K(t,u) = 0 -> solve for u
EX-RK:  A(t,u,ki) = M ki    + K(ti,u0 + dt ∑_{j<i} a_ij * kj) = 0
                  = M ki    + K(ti,ui) = 0 -> solve for ki
where  ui = u0 + dt ∑_{j<i} a_ij * kj for i = 1,…,s
The Jacobian is always M. Compute and store M in EXRungeKuttaStageOperator
"""
function residual!(b::AbstractVector,
  op::EXRungeKuttaStageOperator,
  x::AbstractVector)

  ui = x
  vi = op.vi

  lhs!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. ui = op.u0
  for j = 1:op.i-1
   @. ui = ui  + op.dt * op.a[op.i,j] * op.ki[j]
  end

  rhs = similar(op.u0)
  rhs!(rhs,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. b = b + rhs
  @. b = -1.0 * b
  b
end

function jacobian!(A::AbstractMatrix,
  op::EXRungeKuttaStageOperator,
  x::AbstractVector)
   @. A = op.M
end


function allocate_residual(op::EXRungeKuttaStageOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function allocate_jacobian(op::EXRungeKuttaStageOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end


function update!(op::EXRungeKuttaStageOperator,
  ti::Float64,
  ki::AbstractVector,
  i::Int)
  op.ti = ti
  @. op.ki[i] = ki
  op.i = i
end


function get_mass_matrix!(A::AbstractMatrix,
  odeop::ODEOperator,
  t0::Float64,
  u0::AbstractVector,
  ode_cache)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,odeop,t0,(u0,u0),2,1.0,ode_cache)
  A
end
