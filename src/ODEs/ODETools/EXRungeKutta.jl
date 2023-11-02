"""
Explicit Runge-Kutta ODE solver
"""
struct EXRungeKutta <: ODESolver
  nls::NonlinearSolver
  dt::Float64
  tableau::EXButcherTableau
  function EXRungeKutta(nls::NonlinearSolver, dt, type::Symbol)
    bt = EXButcherTableau(type)
    new(nls, dt, bt)
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
    ki = [similar(u0)]
    M = allocate_jacobian(op,t0,uf,ode_cache)
    get_mass_matrix!(M,op,t0,uf,ode_cache)
    nl_cache = nothing
  else
    ode_cache, vi, ki, nl_cache = cache
  end

  nlop = EXRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,ki,0,a,M)

  i = 1 # for i in 1:s
    # allocate space to store f_i
    if (length(ki) < i)
      push!(ki,similar(u0))
    end

    # solve at stage i
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop,ti,ki,i)
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)

    @. ki[i] = uf
    update!(nlop,ti,ki,i)

  # end

  # update final solution
  tf = t0 + dt
  @. uf = u0 + dt*b[i]*ki[i]
  # for i in 1:s
  # @. uf = uf + dt*b[i]*ki[i]
  # end
  cache = (ode_cache, vi, ki, nl_cache)

  return (uf,tf,cache)


end




mutable struct EXRungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
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
ODE: A(t,u,∂u) = M ∂u/∂t + K(t,u) = 0 -> solve for u
RK:  A(t,u,ki) = M ki    + K(ti,u0 + dt ∑_{j<i} a_ij * kj) = 0 -> solve for ki
               = M ki    + K(ti,ui) = 0
For forward euler, i = 1     -> ui = u0
For other methods, i = 1,…,s -> ui = u0 + dt ∑_{j<i} a_ij * kj
"""
function residual!(b::AbstractVector,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)

  ui = x
  vi = op.vi

  @. vi = 0.0*x
  @. ui = op.u0 # + dt * op.a[op.i,j] * kj
  rhs = similar(op.u0)
  rhs!(rhs,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. ui = 0.0*op.u0
  @. vi = x
  lhs!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. b = b - rhs
  b

end

function jacobian!(A::AbstractMatrix,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  # γ_0^i = 0 (as K is not a function of ki)
  # γ_1^i = 1
  # ui = x
  # vi = op.vi

  # @. ui = op.u0 # this value is irrelevant its jacobian contribution is zero
  # @. vi = x

  # z = zero(eltype(A))
  # fillstored!(A,z)
  # jacobians!(A,op.odeop,op.ti,(ui,vi),(0.0,1.0),op.ode_cache)
  @. A = op.M
end


function allocate_residual(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function allocate_jacobian(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end


function update!(op::EXRungeKuttaStageNonlinearOperator,ti::Float64,ki::AbstractVector,i::Int)
  op.ti = ti
  op.ki = ki
  op.i = i
end



function get_mass_matrix!(A::AbstractMatrix,odeop,t0,u0,ode_cache)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,odeop,t0,(u0,u0),2,1.0,ode_cache)
  A
end
