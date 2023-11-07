using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap
using Test
using Plots

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t) = x -> u(x,t)
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
a(u,v) = ∫(∇(v)⊙∇(u))dΩ
b(v,t) = ∫(v⋅f(t))dΩ
m(ut,v) = ∫(ut⋅v)dΩ

X = TransientMultiFieldFESpace([U,U])
Y = MultiFieldFESpace([V0,V0])

_rhs(t,u,v) = b(v,t) - a(u,v)
_lhs(t,u,v) = ∫( v*u )dΩ

rhs(t,(u1,u2),(v1,v2)) = _rhs(t,u1,v1) + _rhs(t,u2,v2)
lhs(t,(u1,u2),(v1,v2)) = _lhs(t,u1,v1) + _lhs(t,u2,v2)
jac(t,x,(du1,du2),(v1,v2)) = a(du1,v1) + a(du2,v2)
jac_t(t,x,(du1t,du2t),(v1,v2)) = m(du1t,v1) + m(du2t,v2)


op = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,X,Y)

t0 = 0.0
tF = 10.0
dt = 0.001

U0 = U(0.0)
X0 = X(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)
xh0 = interpolate_everywhere([uh0,uh0],X0)

ls = LUSolver()

ode_solver = EXRungeKutta(ls,dt,:EX_FE_1_0_1)

sol_t = solve(ode_solver,op,xh0,t0,tF)

l2(w) = w⋅w

tol = 1.0e-6

errors_rk_fe = []
ts_rk_fe = []
for (xh_tn, tn) in sol_t
  uh_tn = xh_tn[1]
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
  println(el2)

  errors_rk_fe = [errors_rk_fe; el2]
  ts_rk_fe = [ts_rk_fe; tn]
end


plot(ts_rk_fe,errors_rk_fe)
plot!(
  #shape=:auto,
  xlabel="t",
  ylabel="relative error",
  title="RK-FE"
  # yaxis=:log
  )
plot!(show=true)
savefig(string("rk_fe_error_vector"))
