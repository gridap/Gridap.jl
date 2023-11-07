using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap
using Plots
using Test

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
lhs(t,u,v) = ∫( v* (u) )dΩ
rhs(t,u,v) = ∫(v*f(t))dΩ -  ∫(( ∇(v)⊙∇(u) ))dΩ
jac(t,u,du,v) = ∫(( ∇(v)⊙∇(du) ))dΩ
jac_t(t,u,dut,v) = ∫( dut*v )dΩ

op = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V0)

t0 = 0.0
tF = 10.0
dt = 0.001

uh0 = interpolate_everywhere(u(0.0),U(0.0))

ls = LUSolver()
# ode_solver = EXRungeKutta(ls,dt,:EX_FE_1_0_1)
ode_solver = EXRungeKutta(ls,dt,:EX_SSP_3_0_3)

sol_t = solve(ode_solver,op,uh0,t0,tF)


l2(w) = w*w

tol = 1.0e-6

errors_rk_fe = []
ts_rk_fe = []
for (uh_tn, tn) in sol_t
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  println(el2)
  @test el2< tol

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
savefig(string("rk_fe_error_2D"))
