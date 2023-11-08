using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")



using Gridap
using LaTeXStrings
using Plots
using LinearAlgebra
using JLD
using Gridap.Algebra, Gridap.ODEs.ODETools
using GridapSolvers


u(x,t) = x[1]*(1-x[1])*t
# u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t) = x -> u(x,t)
∂tu = ∂t(u)
f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)

n = 2
p = 2
degree = 4*(p+1)
L = 1
dx = L/n
dt = 0.0001
t0 = 0.0
T = 1.0


domain = (0.0, L )
partition = (n)
model  = CartesianDiscreteModel(domain,partition )
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

V = TestFESpace(model,
                ReferenceFE(lagrangian,Float64,p),
                conformity=:H1,
                dirichlet_tags="boundary")
g(x,t::Real) = 0.0
g(t::Real) = x -> g(x,t)
U = TransientTrialFESpace(V,u)



ls = LUSolver()

nls = NewtonRaphsonSolver(ls,1e-10,50)

m(t,u,v) = ∫(v*u)dΩ
lhs(t,u,v) = ∫( v*(u) )dΩ
rhs(t,u,v) = ∫(v*f(t))dΩ - ∫(( ∇(v)⊙∇(u) ))dΩ
jac(t,u,du,v) = ∫(( ∇(v)⊙∇(du) ))dΩ
jac_t(t,u,dut,v) = ∫( dut*v )dΩ



#### SET UP FOR SOLVER
u0 = get_free_dof_values( interpolate_everywhere(u(0),U(0.0)) )
solver = DIRungeKutta(nls,0.001,:BE_1_0_1)
opRK = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
op = Gridap.FESpaces.get_algebraic_operator(opRK)
t0 = 0.0
cache = nothing
uf = copy(u0)


u0 = copy(u_ex)
t0 = tf
###


  # Unpack variables
  dt = solver.dt
  s = solver.tableau.s
  a = solver.tableau.a
  b = solver.tableau.b
  c = solver.tableau.c
  d = solver.tableau.d

  # Create cache if not there
  # if cache === nothing
    ode_cache = Gridap.ODEs.ODETools.allocate_cache(op)
    vi = similar(u0)
    ki = [similar(u0) for i in 1:s]
    nl_cache = nothing
  # else
  #   ode_cache, vi, ki, nl_cache = cache
  # end

  nlop = DIRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,ki,0,a)

  i = 1 #for i in 1:s

    # solve at stage i
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    Gridap.ODEs.ODETools.update!(nlop,ti,ki[i],i)
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)

    Gridap.ODEs.ODETools.update!(nlop,ti,uf,i)

  # end

  # update final solution
  tf = t0 + dt

  @. uf = u0
  for i in 1:s
  @. uf = uf + dt*b[i]*nlop.ki[i]
  end

  cache = (ode_cache, vi, ki, nl_cache)

  return (uf,tf,cache)




u_ex = get_free_dof_values( interpolate_everywhere(u(tf),U(tf)))
println(uf)
println(u_ex)
abs(sum(uf-u_ex))

import Gridap.ODEs.ODETools: RungeKuttaNonlinearOperator
import Gridap.ODEs.ODETools: ODEOperator
mutable struct DIRungeKuttaStageNonlinearOperator <: RungeKuttaNonlinearOperator
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


"""
ODE: A(t,u,∂u) = M ∂u/∂t + K(t,u) = 0 -> solve for u
RK:  A(t,u,ki) = M ki  + K(ti, dt a_ij ki)  + K(ti,u0 + dt ∑_{j<i} a_ij * kj) = 0 -> solve for ki
               = M ki    + K(ti,ui) = 0
For forward euler, i = 1     -> ui = u0
For other methods, i = 1,…,s -> ui = u0 + dt ∑_{j<i} a_ij * kj
"""
function Gridap.Algebra.residual!(b::AbstractVector,op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)

  ui = x
  vi = op.vi

  lhs!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. ui = op.u0
  for j = 1:op.i-1
   @. ui = ui  + op.dt * op.a[op.i,j] * op.ki[j]
  end

  @. ui = ui + op.dt * op.a[op.i,op.i] * x

  rhs = similar(op.u0)
  rhs!(rhs,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. b = b + rhs
  # @. b = -1.0 * b
  println(b)
  b
end

function Gridap.Algebra.jacobian!(A::AbstractMatrix,op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)

  ui = x
  vi = x

  # @. ui = op.u0
  # for j = 1:op.i-1
  #  @. ui = ui  + op.dt * op.a[op.i,j] * op.ki[j]
  # end

  # @. ui = ui + op.dt * op.a[op.i,op.i] * x

  # @. vi = x


  z = zero(eltype(A))
  LinearAlgebra.fillstored!(A,z)
  Gridap.ODEs.ODETools.jacobians!(A,op.odeop,op.ti,(ui,vi),(op.dt*op.a[op.i,op.i],1.0),op.ode_cache)

end


function Gridap.Algebra.allocate_residual(op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function Gridap.Algebra.allocate_jacobian(op::DIRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end


function Gridap.ODEs.ODETools.update!(op::DIRungeKuttaStageNonlinearOperator,ti::Float64,ki::AbstractVector,i::Int)
  op.ti = ti
  @. op.ki[i] = ki
  op.i = i
end
