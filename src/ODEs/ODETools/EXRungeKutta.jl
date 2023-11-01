"""
Explicit Runge-Kutta ODE solver.

This struct defines an ODE solver for the system of ODEs

    M(u,t)du/dt = f(u,t)

where `f` is a nonlinear function of `u` and `t` that will be treated explicitly.
  The ODE is solved using an explicit Runge-Kutta method.
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
    ki = [similar(u0)]
    nl_cache = nothing
  else
    ode_cache, vi, ki, nl_cache = cache
  end

  nlop = EXRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,ki,0,a)

  for i in 1:s
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



  end

  # update
  @. uf = u0
  for i in 1:s
  @. uf = uf + dt*b[i]*ki[i]
  end
  cache = (ode_cache, vi, ki, nl_cache)
  tf = t0 + dt

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
end



mutable struct EXRungeKuttaUpdateNonlinearOperator <: RungeKuttaNonlinearOperator
  odeop::ODEOperator
  ti::Float64
  dt::Float64
  u0::AbstractVector
  ode_cache
  vi::AbstractVector
end

function residual!(b::AbstractVector,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  # A(t,ui,∂ui/∂t) = ∂ui/∂t - ∑_{j<i} a_ij * f(tj,uj) = 0
  # b = [∂ui/∂t ]
  # b - ∑_{j<i} a_ij * f(tj,uj) = 0

  # lhs!(b,op,x)
  # rhs!(op,x)
  # @. b = b - op.rhs
  # b
  vi = op.vi
  @. vi = (x-op.u0)/(op.dt)
  ui = op.u0 #op.a[op.i,op.i] * x
  residual!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)

end

function jacobian!(A::AbstractMatrix,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  # A(t,ui,∂ui/∂t) = ∂ui/∂t - ∑_{j<i} a_ij * f(tj,uj) = 0
  # γ_0^i = 0 (as A is not a function of u_i)
  # γ_1^i = 1/δt
  ui = x
  vi = op.vi
  @. vi = (x-op.u0)/(op.dt)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.ti,(ui,vi),(0.0,1.0),op.ode_cache)
end


function allocate_residual(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function allocate_jacobian(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end


# function rhs!(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
#   v = op.vi
#   @. v = (x-op.u0)/(op.dt)
#   u = op.a[op.i,op.i] * x # == zero for explicit
#   # if (op.i>1)
#   #   @. u += op.ui
#   # end
#   rhs!(op.rhs,op.odeop,op.ti,(u,v),op.ode_cache)

# end

# function lhs!(b::AbstractVector,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
#   ui = x
#   vi = op.vi
#   @. vi = (x-op.u0)/(op.dt)
#   lhs!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)
# end

function update!(op::EXRungeKuttaStageNonlinearOperator,ti::Float64,ki::AbstractVector,i::Int)
  op.ti = ti
  op.ki = ki
  op.i = i
end
