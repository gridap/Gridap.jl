using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap
using LaTeXStrings
using Plots
using LinearAlgebra
using JLD
using Gridap.Algebra, Gridap.ODEs.ODETools



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

m(t,u,v) = ∫(v*u)dΩ
lhs(t,u,v) = ∫( v*(u) )dΩ
rhs(t,u,v) = ∫(v*f(t))dΩ - ∫(( ∇(v)⊙∇(u) ))dΩ
# res(t,u,v) = lhs(t,u,v) - rhs(t,u,v)
jac(t,u,du,v) = ∫(( ∇(v)⊙∇(du) ))dΩ
jac_t(t,u,dut,v) = ∫( dut*v )dΩ



#### SET UP FOR SOLVER
u0 = get_free_dof_values( interpolate_everywhere(u(0),U(0.0)) )
# solver = EXRungeKutta(ls,0.001,:EX_FE_1_0_1)
solver = EXRungeKutta(ls,0.001,:EX_SSP_3_0_3)
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
    ode_cache = Gridap.ODEs.ODETools.allocate_cache(op)
    vi = similar(u0)
    ki = [similar(u0)]
    M = allocate_jacobian(op,t0,uf,ode_cache)
    get_mass_matrix!(M,op,t0,uf,ode_cache)
    nl_cache = nothing

  nlop = EXRungeKuttaStageNonlinearOperator(op,t0,dt,u0,ode_cache,vi,ki,0,a,M)
  println(nlop.ti)
  println(nlop.u0)
  println(nlop.ki)
  println(nlop.M)

  for i in 1:s
    # allocate space to store k_i
    if (length(ki) < i)
      push!(ki,similar(u0))
    end

    # solve at stage i
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    Gridap.ODEs.ODETools.update!(nlop,ti,ki,i)
    nl_cache = solve!(uf,solver.nls,nlop,nl_cache)

    @. ki[i] = uf
    Gridap.ODEs.ODETools.update!(nlop,ti,ki,i)


  end

  # update
  @. uf = u0
  for i in 1:s
  @. uf = uf + dt*b[i]*ki[i]
  end
  cache = nothing #(ode_cache, vi, ki, nl_cache)
  tf = t0 + dt

  # return (uf,tf,cache)


u_ex = get_free_dof_values( interpolate_everywhere(u(tf),U(tf)))
println(uf)
println(u_ex)
abs(sum(uf-u_ex))

import Gridap.ODEs.ODETools: RungeKuttaNonlinearOperator
import Gridap.ODEs.ODETools: ODEOperator
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
function Gridap.Algebra.residual!(b::AbstractVector,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)

  ui = x
  vi = op.vi

  lhs!(b,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. ui = op.u0
  for j = 1:op.i-1
   @. ui = ui  + dt * op.a[op.i,j] * op.ki[j]
  end

  rhs = similar(op.u0)
  rhs!(rhs,op.odeop,op.ti,(ui,vi),op.ode_cache)

  @. b = b + rhs
  @. b = -1.0 * b
  b

end




function Gridap.Algebra.jacobian!(A::AbstractMatrix,op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  # γ_0^i = 0 (as K is not a function of ki)
  # γ_1^i = 1

  @. A = op.M

  # ui = x
  # vi = op.vi
  # @. ui = op.u0 # this value is irrelevant its jacobian contribution is zero
  # @. vi = x
  # z = zero(eltype(A))
  # LinearAlgebra.fillstored!(A,z)
  # Gridap.ODEs.ODETools.jacobians!(A,op.odeop,op.ti,(ui,vi),(0.0,1.0),op.ode_cache)

end


function Gridap.Algebra.allocate_residual(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_residual(op.odeop,op.ti,x,op.ode_cache)
end

function Gridap.Algebra.allocate_jacobian(op::EXRungeKuttaStageNonlinearOperator,x::AbstractVector)
  allocate_jacobian(op.odeop,op.ti,x,op.ode_cache)
end


function Gridap.ODEs.ODETools.update!(op::EXRungeKuttaStageNonlinearOperator,ti::Float64,ki::AbstractVector,i::Int)
  op.ti = ti
  op.ki = ki
  op.i = i
end


function get_mass_matrix!(A::AbstractMatrix,odeop,t0,u0,ode_cache)
  z = zero(eltype(A))
  LinearAlgebra.fillstored!(A,z)
  Gridap.ODEs.ODETools.jacobian!(A,odeop,t0,(u0,u0),2,1.0,ode_cache)
  A
end

# M = allocate_jacobian(op,t0,uf,ode_cache) #allocate_jacobian(nlop.odeop,nlop.ti,uf,nlop.ode_cache)
# get_mass_matrix!(M,op,t0,uf,ode_cache)


# ki[i] = Vector{Float64}([3, 6, 10])
# b = similar(nlop.u0)
# mul!(b,nlop.M,ki[i])

# vi = ones(size(nlop.u0))
# ui = ki[i] # + dt * op.a[op.i,j] * kj
# _lhs = similar(nlop.u0)
# Gridap.ODEs.ODETools.lhs!(_lhs,nlop.odeop,nlop.ti,(ui,vi),nlop.ode_cache)
# print( sum( b - _lhs)  )

# ui = 100*ones(size(nlop.u0))
# vi = nlop.vi
# frhs = similar(nlop.u0)
# Gridap.ODEs.ODETools.rhs!(frhs,nlop.odeop,nlop.ti,(ui,vi),nlop.ode_cache)

# b = similar(nlop.u0)
# Gridap.ODEs.ODETools.lhs!(b,nlop.odeop,nlop.ti,(ui,vi),nlop.ode_cache)
