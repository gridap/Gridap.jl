module ConstantFEOperatorsTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.ODEs.ODETools
using Gridap.ODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator
using LineSearches: BackTracking

import Gridap: ∇
import Gridap.ODEs.TransientFETools: ∂t

θ = 1.0

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]
u(t::Real) = x -> u(x,t)
∂tu = ∂t(u)

f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 2

reffe = ReferenceFE(lagrangian,Float64,order)
V0 = FESpace(
  model,
  reffe,
  conformity=:H1,
  dirichlet_tags="boundary"
)
U = TransientTrialFESpace(V0,u)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

#
a(u,v) = ∫(∇(v)⋅∇(u))dΩ
b(v) = ∫(v*f(0.0))dΩ
m(ut,v) = ∫(ut*v)dΩ

op = TransientConstantFEOperator(m,a,b,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ls = LUSolver()
ode_solver = ThetaMethod(ls,dt,θ)
sol_t = solve(ode_solver,op,uh0,t0,tF)

l2(w) = w*w

tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

end #module
