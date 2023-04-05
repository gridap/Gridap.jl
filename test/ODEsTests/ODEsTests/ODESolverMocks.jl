using Gridap.Algebra: residual
using Gridap.Algebra: jacobian
import Gridap.Algebra: NonlinearSolver
import Gridap.Algebra: NonlinearOperator
import Gridap.Algebra: solve!
import Gridap.ODEs.ODETools: solve_step!
import Gridap.ODEs.ODETools: ODESolver
import Gridap.ODEs.ODETools: zero_initial_guess
import Gridap.ODEs.ODETools: residual!
import Gridap.ODEs.ODETools: jacobian!
import Gridap.ODEs.ODETools: solve!
import Gridap.ODEs.ODETools: allocate_residual
import Gridap.ODEs.ODETools: allocate_jacobian

struct OperatorMock <: NonlinearOperator
  odeop
  tf::Float64
  dt::Float64
  u0::AbstractVector
  cache
end

function OperatorMock(odeop::ODEOperator,tf::Real,dt::Real,u0::AbstractVector)
  cache = nothing
  OperatorMock(odeop,tf,dt,u0,cache)
end

function residual!(b::AbstractVector,op::OperatorMock,x::AbstractVector)
  uf = x
  uf_t = (x-op.u0)/op.dt
  residual!(b,op.odeop,op.tf,(uf,uf_t),op.cache)
end

function jacobian!(A::AbstractMatrix,op::OperatorMock,x::AbstractVector)
  uf = x
  uf_t = (x-op.u0)/op.dt
  fill!(A,0.0)
  # jacobian!(A,op.odeop,op.tf,(uf,uf_t),1,1.0,op.cache)
  # jacobians!(A,op.odeop,op.tf,(uf,uf_t),(1/op.dt),op.cache)
  jacobians!(A,op.odeop,op.tf,(uf,uf_t),(1.0,1.0/op.dt),op.cache)
end

function allocate_residual(op::OperatorMock,x::AbstractVector)
  allocate_residual(op.odeop,op.tf,x,op.cache)
end

function allocate_jacobian(op::OperatorMock,x::AbstractVector)
  allocate_jacobian(op.odeop,op.tf,x,op.cache)
end

function zero_initial_guess(op::OperatorMock)
  x0 = similar(op.u0)
  fill!(x0,zero(eltype(x0)))
  x0
end

struct NLSolverMock <: NonlinearSolver
end

function solve!(x::AbstractVector,nls::NLSolverMock,nlop::NonlinearOperator,cache::Nothing)
  r = residual(nlop,x)
  J = jacobian(nlop,x)
  dx = inv(Matrix(J))*(-r)
  x.= x.+dx
  cache = (r,J,dx)
end

function solve!(x::AbstractVector,nls::NLSolverMock,nlop::NonlinearOperator,cache)
  r, J, dx = cache
  residual!(r, nlop, x)
  jacobian!(J, nlop, x)
  dx = inv(Matrix(J))*(-r)
  x.= x.+dx
end

struct ODESolverMock <: ODESolver
  nls::NLSolverMock
  dt::Float64
end

function solve_step!(
  uf::AbstractVector,solver::ODESolverMock,op::ODEOperator,u0::AbstractVector,t0::Real, cache) # -> (uF,tF)

  dt = solver.dt
  tf = t0+dt
  if (cache == nothing)
    ode_cache = allocate_cache(op)
  else
    ode_cache, nl_cache = cache
    update_cache!(ode_cache,op,tf)
  end

  nlop = OperatorMock(op,tf,dt,u0,ode_cache)

  if (cache==nothing)
    nl_cache = solve!(uf,solver.nls,nlop)
  else
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)
  end

  cache = ode_cache, nl_cache

  return (uf, tf, cache)
end
