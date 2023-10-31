using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap
# using LaTeXStrings
using Plots
# using JLD


function l2(u,Ω,p)
  dΩ_error = Measure(Ω, 10*(p + 1))
  return sqrt(sum( ∫(u ⊙ u)*dΩ_error))
end


u(x,t) = x[1]*(1-x[1])*t
u(t) = x -> u(x,t)
∂tu = ∂t(u)
f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)


n = 4
p = 2
degree = 4*(p+1)
L = 1
dx = L/n
dt = 0.0001
t0 = 0.0
T = 1.0


domain = (0.0, L)
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
U = TransientTrialFESpace(V,g)

u0 = interpolate_everywhere(u(0),U(0.0))


ls = LUSolver()

a(t,u,v) = ∫(( ∇(v)⊙∇(u) ))dΩ
m(u,v) = ∫(v*u)dΩ
b(v,t) = ∫(v*f(t))dΩ
lhs(t,u,v) = m(u,v)
rhs(t,u,v) = b(v,t) - a(t,u,v)
res(t,u,v) = a(t,u,v) + m(∂t(u),v) - b(v,t)
jac(t,u,du,v) = a(t,du,v)
jac_t(t,u,dut,v) = m(dut,v)



#### Solve with standard BE
be = BackwardEuler(ls,0.001)
op_be = TransientFEOperator(res,jac,jac_t,U,V)
sol_be = solve(be,op_be,u0,t0,T)

errors = []
ts = []
for (uh,t) in sol_be
  u_ex = interpolate_everywhere(u(t),U(t))
  e = l2(uh-u_ex,Ω,p)/l2(u_ex,Ω,p) # relative error
  errors = [errors; e]
  ts = [ts; t]
end

plot()
plot!(ts,errors)
plot!(
  #shape=:auto,
  xlabel="t",
  ylabel="relative error",
  title="BackwardEuler"
  # yaxis=:log
  )
plot!(show=true)
savefig(string("be_error"))


#### Solve with standard RungeKutta with BE table
rk = RungeKutta(ls,ls,0.001,:BE_1_0_1)
opRK = TransientFEOperator(res,jac,jac_t,U,V) #TransientRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
sol_rk = solve(rk,opRK,u0,t0,T)


errors_rk = []
ts_rk = []
for (uh,t) in sol_rk
  u_ex = interpolate_everywhere(u(t),U(t))
  e = l2(uh-u_ex,Ω,p)/l2(u_ex,Ω,p) # relative error
  errors_rk = [errors_rk; e]
  ts_rk = [ts_rk; t]
end

plot()
plot!(ts_rk,errors_rk)
plot!(
  #shape=:auto,
  xlabel="t",
  ylabel="relative error",
  title="RK-BE"
  # yaxis=:log
  )
plot!(show=true)
savefig(string("rk_be"))




#### Solve with standard FE
fe = ForwardEuler(ls,0.001)
op_fe = TransientFEOperator(res,jac,jac_t,U,V)
sol_fe = solve(fe,op_fe,u0,t0,T)

errors_fe = []
ts_fe = []
for (uh,t) in sol_fe
  u_ex = interpolate_everywhere(u(t),U(t))
  e = l2(uh-u_ex,Ω,p)/l2(u_ex,Ω,p) # relative error
  errors_fe = [errors_fe; e]
  ts_fe = [ts_fe; t]
end

plot()
plot!(ts_fe,errors_fe)
plot!(
  #shape=:auto,
  xlabel="t",
  ylabel="relative error",
  title="ForwardEuler"
  # yaxis=:log
  )
plot!(show=true)
savefig(string("fe_error"))

#### Solve with standard EXRungeKutta with FE table
rk_fe = EXRungeKutta(ls,0.001,:EX_FE_1_0_1)
opRK_fe = TransientFEOperator(res,jac,jac_t,U,V) #TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
sol_rk_fe = solve(rk_fe,opRK_fe,u0,t0,T)

errors_rk_fe = []
ts_rk_fe = []
for (uh,t) in sol_rk_fe
  u_ex = interpolate_everywhere(u(t),U(t))
  e = l2(uh-u_ex,Ω,p)/l2(u_ex,Ω,p) # relative error
  errors_rk_fe = [errors_rk_fe; e]
  ts_rk_fe = [ts_rk_fe; t]
end

plot()
plot!(ts_rk_fe,errors_rk_fe)
plot!(
  #shape=:auto,
  xlabel="t",
  ylabel="relative error",
  title="RK-FE"
  # yaxis=:log
  )
plot!(show=true)
savefig(string("rk_fe"))