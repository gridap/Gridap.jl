module VectorHeatEquationTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.FESpaces: get_algebraic_operator

θ = 1.0

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t::Real) = x -> u(x,t)
v(x) = t -> u(x,t)
f(t) = x -> ∂t(u)(t)(x)-Δ(u(t))(x)

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
a(u,v) = ∫(∇(v)⊙∇(u))dΩ
b(v,t) = ∫(v⋅f(t))dΩ
m(ut,v) = ∫(ut⋅v)dΩ

X = TransientMultiFieldFESpace([U,U])
Y = MultiFieldFESpace([V0,V0])

_res(t,u,v) = a(u,v) + m(∂t(u),v) - b(v,t)

res(t,(u1,u2),(v1,v2)) = _res(t,u1,v1) + _res(t,u2,v2)
jac(t,x,(du1,du2),(v1,v2)) = a(du1,v1) + a(du2,v2)
jac_t(t,x,(du1t,du2t),(v1,v2)) = m(du1t,v1) + m(du2t,v2)

op = TransientFEOperator(res,jac,jac_t,X,Y)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
X0 = X(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)
xh0 = interpolate_everywhere([uh0,uh0],X0)

ls = LUSolver()
# using Gridap.Algebra: NewtonRaphsonSolver
# nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
ode_solver = ThetaMethod(ls,dt,θ)

sol_t = solve(ode_solver,op,xh0,t0,tF)

l2(w) = w⋅w


tol = 1.0e-6
_t_n = t0

result = Base.iterate(sol_t)

for (xh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  uh_tn = xh_tn[1]
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

#copy test
all_sol = [ (copy(xh_tn[1]), tn) for (xh_tn, tn) in sol_t ]
all_el2 = [ sqrt(sum( ∫(l2( u(tn) - uhc_tn ))dΩ )) for (uhc_tn,tn) in all_sol ]
@test all( all_el2 .< tol )

end #module
