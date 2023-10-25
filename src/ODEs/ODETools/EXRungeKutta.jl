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
    fi = [similar(u0)]
    nl_cache = nothing
  else
    ode_cache, vi, fi, nl_cache = cache
  end

  nlop = EXRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,fi,0,a)

  for i in 1:s
     # allocate space to store the RHS at i
    if (length(fi) < i)
      push!(fi,similar(u0))
    end

    # solve at stage i
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop,ti,fi,i)
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)
    fi[i] = get_fi(uf,nlop,nl_cache)

  end

  # update
  uf = u0
  for i in 1:s
    uf = uf + dt*b[i]*fi[i]
  end
  cache = (ode_cache, vi, fi, nl_cache)
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
  fi::AbstractVector
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
  # A(t,ui,∂ui/∂t) = ∂ui/∂t - ∑_{j<i} a_ij * f(uj,tj) = 0
  # b = [∂u/∂t ]
  # Res_ij = - a_ij * f(uj,ti)
  # b + ∑_{j<i} Res_ij = 0
  ui = x
  vi = op.vi
  vi = (x-op.u0)/(op.dt)
  residual!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)
  for j in 1:op.i-1
    b .= b - op.a[op.i,j]* op.fi[j]
  end
  b
end
# function residual!(b::AbstractVector,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
#   vi = op.vi
#   @. vi = (x-op.u0)/op.dt
#   residual!(b,op.odeop,op.ti,(op.u0,vi),op.ode_cache) # in FE, use u0 not ui
# end


function jacobian!(A::AbstractMatrix,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  # @assert (abs(op.a[op.i,op.i]) > 0.0)
  ui = x
  vi = op.vi
  vi = (x-op.u0)/(op.dt)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.ti,(ui,vi),(0.0,1.0),op.ode_cache)
end
# function jacobian!(A::AbstractMatrix,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
#   # @assert (abs(op.a[op.i,op.i]) > 0.0)
#   # ui = x # this line not in FE
#   vi = op.vi
#   @. vi = (x-op.u0)/(op.dt)
#   z = zero(eltype(A))
#   fillstored!(A,z)
#   jacobians!(A,op.odeop,op.ti,(op.u0,vi),(0,1.0/op.dt),op.ode_cache) # in FE, use u0 not ui
# end

function allocate_residual(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function allocate_jacobian(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end


function get_fi(x::AbstractVector, op::EXRungeKuttaStageNonlinearOperator, cache::Nothing)
  ui = x
  vi = op.vi
  vi=zero(x)
  b=similar(x)
  residual!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)
  (vi-b) # store fi for future stages
end
function get_fi(x::AbstractVector, op::EXRungeKuttaStageNonlinearOperator, cache)
  ui = x
  vi = op.vi
  vi=zero(x)
  residual!(cache.b,op.odeop,op.ti,(ui,vi),op.ode_cache)
  (vi-cache.b) # store fi for future stages
end


function update!(op::EXRungeKuttaStageNonlinearOperator,ti::Float64,fi::AbstractVector,i::Int)
  op.ti = ti
  op.fi = fi
  op.i = i
end




"""
EXRungeKuttaStageNonlinearOperator <: NonlinearOperator

Nonlinear operator that represents the Runge-Kutta nonlinear operator at a
given time step and stage, i.e., A(t,u_i,(u_i-u_n)/dt)
"""
# mutable struct EXRungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
#   odeop::ODEOperator
#   ti::Float64
#   dt::Float64
#   u0::AbstractVector
#   ode_cache
#   vi::AbstractVector
#   ui::Vector{AbstractVector}
#   rhs::AbstractVector
#   i::Int
#   a::Matrix{Float64}
# end

# """
# EXRungeKuttaUpdateNonlinearOperator <: NonlinearOperator

# Nonlinear operator that represents the Runge-Kutta nonlinear operator at the
# final updated of a given time step, i.e., A(t,u_f,(u_f-u_n)/dt)
# """
# mutable struct EXRungeKuttaUpdateNonlinearOperator <: RungeKuttaNonlinearOperator
#   odeop::ODEOperator
#   ti::Float64
#   dt::Float64
#   u0::AbstractVector
#   ode_cache
#   vi::AbstractVector
#   ui::Vector{AbstractVector}
#   rhs::AbstractVector
#   s::Int
#   b::Vector{Float64}
# end

# """
# residual!(b,op::EXRungeKuttaStageNonlinearOperator,x)

# Compute the residual of the Runge-Kutta nonlinear operator at stage i.
# ```math
# A(t,ui,∂ui/∂t) = ∂ui/∂t - a_ii * f(ui,ti) - ∑_{j<i} a_ij * f(uj,tj) = 0
# ```

# Uses the vector b as auxiliar variable to store the residual of the left-hand side of
# the i-th stage ODE operator, then adds the corresponding contribution from right-hand side
# at all earlier stages.
# ```math
# b = M(ui,ti)∂u/∂t
# b - f(∑_{j<=i} a_ij * uj,tj) = 0
# ```
# """
# function residual!(b::AbstractVector,
#   op::EXRungeKuttaStageNonlinearOperator,
#   x::AbstractVector)
#   rhs!(op,x) ### this is different to below
#   lhs!(b,op,x)
#   @. b = b - op.rhs
#   b
# end

# """
# residual!(b,op::EXRungeKuttaUpdateNonlinearOperator,x)

# Compute the residual of the final update of the Runge-Kutta nonlinear operator.
# ```math
# A(t,uf,∂uf/∂t) = ∂uf/∂t - ∑_{i<=s} b_i * f(ui,ti) = 0
# ```

# Uses the vector b as auxiliar variable to store the residual of the update ODE
# operator (e.g. identity or mass matrix), then adds the corresponding contribution from earlier stages.
# ```math
# b = M(uf,tf)[∂u/∂t]
# b - f(∑_{i<=s} bi * ui,ti) = 0
# ```
# """
# function residual!(b::AbstractVector,
#   op::EXRungeKuttaUpdateNonlinearOperator,
#   x::AbstractVector)
#   lhs!(b,op,x)
#   rhs!(op,x)
#   @. b = b - op.rhs
#   b
# end

# function jacobian!(A::AbstractMatrix,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
#   # @assert (abs(op.a[op.i,op.i]) > 0.0)
#   ui = x
#   vi = op.vi
#   @. vi = (x-op.u0)/(op.dt)
#   z = zero(eltype(A))
#   fillstored!(A,z)
#   jacobians!(A,op.odeop,op.ti,(ui,vi),(op.a[op.i,op.i],1.0/op.dt),op.ode_cache)
#   # do not need to compute J0 as a_ii = 0
# end

# function jacobian!(A::AbstractMatrix,op::EXRungeKuttaUpdateNonlinearOperator,x::AbstractVector)
#   uf = x
#   vf = op.vi
#   @. vf = (x-op.u0)/(op.dt)
#   z = zero(eltype(A))
#   fillstored!(A,z)
#   jacobian!(A,op.odeop,op.ti,(uf,vf),1.0,1.0/(op.dt),op.ode_cache)
#   # I have changed this input from hard coded 2 -> 1
# end

# function jacobian!(A::AbstractMatrix,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector,
#   i::Integer,γᵢ::Real)
#   ui = x
#   vi = op.vi
#   @. vi = (x-op.u0)/(op.dt)
#   z = zero(eltype(A))
#   fillstored!(A,z)
#   jacobian!(A,op.odeop,op.ti,(ui,vi),i,γᵢ,op.ode_cache)
# end

# function allocate_residual(op::RungeKuttaNonlinearOperator,x::AbstractVector)
#   allocate_residual(op.odeop,op.ti,x,op.ode_cache)
# end

# function allocate_jacobian(op::RungeKuttaNonlinearOperator,x::AbstractVector)
#   allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
# end

# function zero_initial_guess(op::RungeKuttaNonlinearOperator)
#   x0 = similar(op.u0)
#   fill!(x0,zero(eltype(x0)))
#   x0
# end

# function rhs!(op::EXRungeKuttaStageNonlinearOperator, x::AbstractVector)
#   u = x
#   u .= 0.0 # changed this based on Santi comment on Oriol code
#   for j in 1:op.i-1 # changed to i-1 for explicit methods
#     @. u = u + op.a[op.i,j] * op.ui[j]
#   end
#   v = op.vi
#   @. v = (x-op.u0)/(op.dt)
#   rhs!(op.rhs,op.odeop,op.ti,(u,v),op.ode_cache)
# end

# function rhs!(op::EXRungeKuttaUpdateNonlinearOperator, x::AbstractVector)
#   u = x
#   u .= 0.0 # changed this based on Santi comment on Oriol code
#   for i in 1:op.s
#     @. u = u + op.b[i] * op.ui[i]
#   end
#   v = op.vi
#   @. v = (x-op.u0)/(op.dt)
#   rhs!(op.rhs,op.odeop,op.ti,(u,v),op.ode_cache)
# end

# function lhs!(b::AbstractVector, op::RungeKuttaNonlinearOperator, x::AbstractVector)
#   u = x
#   v = op.vi
#   @. v = (x-op.u0)/(op.dt)
#   lhs!(b,op.odeop,op.ti,(u,v),op.ode_cache)
# end

# function update!(op::RungeKuttaNonlinearOperator,ti::Float64,ui::AbstractVector,i::Int)
#   op.ti = ti
#   op.ui = ui
#   op.i = i
# end

# function _mass_matrix!(A,odeop,t,ode_cache,u)
#   z = zero(eltype(A))
#   fillstored!(A,z)
#   jacobian!(A,odeop,t,(u,u),2,1.0,ode_cache)
# end
