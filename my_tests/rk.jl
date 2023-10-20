using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")
Pkg.add("Plots")


using Gridap
using LaTeXStrings
using Plots
using JLD


n = 16
p = 1
degree = 4*(p+1)
L = 1
dx = L/n
dt = 0.0001
nT = 100
t0 = 0.0
T = 1.0


domain = (0.0, L, 0.0, L)
partition = (n, n)
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

u0 = interpolate_everywhere(0.0,U(0.0))


ls = LUSolver()
f(t) = sin(2*π*t)
κ(t) = 1.0 + 0.95*sin(2*π*t) #1.0 + 0.01*sin(2*π*t)

a(t,u,v) = ∫(κ(t) *( ∇(v)⊙∇(u) ))dΩ
m(u,v) = ∫(v*u)dΩ
b(v,t) = ∫(v*f(t))dΩ
lhs(t,u,v) = m(u,v)
rhs(t,u,v) = b(v,t) - a(t,u,v)
res(t,u,v) = a(t,u,v) + m(∂t(u),v) - b(v,t)
jac(t,u,du,v) = a(t,du,v)
jac_t(t,u,dut,v) = m(dut,v)

u_fe = []
u_rk_ex = []
ts = []


#### Solve with EXRungeKutta
rk_ex = EXRungeKutta(ls,ls,0.001,:EX_FE_1_0_1)
opRK_ex = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
opRK_ex = TransientEXRungeKuttaFEOperator(m,rhs,U,V)
sol_rk_ex = solve(rk_ex,opRK_ex,u0,t0,T)




createpvd("my_tests/transient_sol/rk_ex") do pvd
  for (uh,t) in sol_rk_ex
    if ( mod(t,0.1)  < 1e-4  )
      # u_rk_ex = [u_rk_ex; uh]
      # ts = [ts; t]

      pvd[t] = createvtk(Ω,"my_tests/transient_sol/rk_ex_$t"*".vtu",cellfields=["u"=>uh])
      println(t)
    end
  end
end


#### Solve with standard RungeKutta with FE table
rk = RungeKutta(ls,ls,0.001,:BE_1_0_1)
opRK = TransientRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
# opRK = TransientRungeKuttaFEOperator(m,rhs,U,V)
sol_rk = solve(rk,opRK,u0,t0,T*5)

createpvd("my_tests/transient_sol/rk_og") do pvd
  for (uh,t) in sol_rk
    # if ( mod(t,0.1)  < 1e-4  )
      pvd[t] = createvtk(Ω,"my_tests/transient_sol/rk_og_$t"*".vtu",cellfields=["u"=>uh])
      println(t)
    # end
  end
end

#### Solve with forward Euler method in Gridap
fe = ForwardEuler(ls,dt) # ThetaMethod(ls,dt,0.0)
op_fe = TransientFEOperator(res,jac,jac_t,U,V)
sol_fe = solve(fe,op_fe,u0,t0,T)


createpvd("my_tests/transient_sol/rk_fe") do pvd
  for (uh,t) in sol_fe
    if ( mod(t,0.1)  < 1e-4  )
      # u_fe = [u_fe; uh]

      pvd[t] = createvtk(Ω,"my_tests/transient_sol/rk_fe_$t"*".vtu",cellfields=["u"=>uh])
      println(t)
    end
  end
end


function l2(u,Ω,p)
  dΩ_error = Measure(Ω, 10*(p + 1))
  return sqrt(sum( ∫(u ⊙ u)*dΩ_error))
end


errors = zeros(size(u_rk_ex))
for i in 1:length(u_rk_ex)
  errors[i] = l2( (u_rk_ex[i]-u_fe[i] ) ,Ω,p)
end

plot()
plot!(ts,errors,lw=2)
plot!(
  #shape=:auto,
  xlabel=L"$t$",
  ylabel="error",
  # yaxis=:log,
  ylimits=(0.0115,0.0125),
  )
plot!(show=true)
savefig(string("l2_error_p",p))



###############################################################
# manufactured soln




using Gridap
using LaTeXStrings
using Plots
using JLD

import Gridap.FESpaces: get_algebraic_operator


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


ode_op = get_algebraic_operator(op_be)
get_matrix(ode_op)
ode_op.matrix


createpvd("transient_sol/be") do pvd
  for (uh,t) in sol_be
    # if ( mod(t,0.1)  < 1e-4  )
      pvd[t] = createvtk(Ω,"transient_sol/be_$t"*".vtu",cellfields=["u"=>uh])
      println(t)
    # end
  end
end

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
opRK = TransientRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
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



createpvd("transient_sol/rk_be") do pvd
  for (uh,t) in sol_rk
    # if ( mod(t,0.1)  < 1e-4  )
      pvd[t] = createvtk(Ω,"transient_sol/rk_be_$t"*".vtu",cellfields=["u"=>uh])
      println(t)
    # end
  end
end
