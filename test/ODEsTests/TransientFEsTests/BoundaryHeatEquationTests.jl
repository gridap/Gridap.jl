module BoundaryHeatEquationTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.FESpaces: get_algebraic_operator

import Gridap: ∇
import Gridap.ODEs.TransientFETools: ∂t

θ = 0.2

# Analytical functions
# u(x,t) = (x[1]+x[2])*t
# u(x,t) = (2*x[1]+x[2])*t
u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t::Real) = x -> u(x,t)
v(x) = t -> u(x,t)
∂tu(t) = x -> ForwardDiff.derivative(v(x),t)
∂tu(x,t) = ∂tu(t)(x)
∂t(::typeof(u)) = ∂tu
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
  dirichlet_tags=[1,2,3,4,5,6]
)
U = TransientTrialFESpace(V0,u)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

neumanntags = [7,8]
Γ = BoundaryTriangulation(model,tags=neumanntags)
dΓ = Measure(Γ,degree)
nb = get_normal_vector(Γ)

#
a(u,v) = ∫(∇(v)⋅∇(u))dΩ
b(v,t) = ∫(v*f(t))dΩ
m(ut,v) = ∫(ut*v)dΩ
b_Γ(v,t) = ∫(v*(∇(u(t))⋅nb))dΓ

res(t,u,v) = a(u,v) + m(∂t(u),v) - b(v,t) - b_Γ(v,t)
jac(t,u,du,v) = a(du,v)
jac_t(t,u,dut,v) = m(dut,v)

op = TransientFEOperator(res,jac,jac_t,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ls = LUSolver()
using Gridap.Algebra: NewtonRaphsonSolver
# nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
ode_solver = ThetaMethod(ls,dt,θ)

sol_t = solve(ode_solver,op,uh0,t0,tF)

# Juno.@enter Base.iterate(sol_t)

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

all_sol = [ (copy(uh_tn), tn) for (uh_tn, tn) in sol_t ]
all_el2 = [ sqrt(sum( ∫(l2( u(tn) - uhc_tn ))dΩ )) for (uhc_tn,tn) in all_sol ]
@test all( all_el2 .< tol )

end #module
