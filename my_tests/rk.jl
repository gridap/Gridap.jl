using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")


using Gridap
using LaTeXStrings
using Plots
using JLD


n = 64
p = 1
degree = 4*(p+1)
L = 1
dx = L/n
dt = 0.05
nT = 200
t0 = 0.0
T = dt*nT


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

f(t) = sin(π*t)
κ(t) = 1.0 + 0.95*sin(2π*t)


a(t,u,v) = ∫(κ(t) *( ∇(v)⊙∇(u) ))dΩ
m(u,v) = ∫(v*u)dΩ
b(v,t) = ∫(v*f(t))dΩ
lhs(t,u,v) = m(u,v)
rhs(t,u,v) = b(v,t) - a(t,u,v)
jac(t,u,du,v) = a(t,du,v)
jac_t(t,u,dut,v) = m(dut,v)


rk_ex = EXRungeKutta(ls,ls,0.005,:EX_FE_1_0_1)
opRK_ex = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
sol_rk_ex = solve(rk_ex,opRK_ex,u0,t0,T)




for (uh,t) in sol_rk_ex
  println(uh)
end


createpvd("my_tests/transient_sol/rk") do pvd
  for (uh,t) in sol_rk_ex
    pvd[t] = createvtk(Ω,"my_tests/transient_sol/rk_$t"*".vtu",cellfields=["u"=>uh])
  end
end




a(t,u,v) = ∫(κ(t) *( ∇(v)⊙∇(u) ))dΩ
m(u,v) = ∫(v*u)dΩ
b(v,t) = ∫(v*f(t))dΩ
lhs(t,u,v) = m(u,v)
rhs(t,u,v) = b(v,t) - a(t,u,v)
jac(t,u,du,v) = a(t,du,v)
jac_t(t,u,dut,v) = m(dut,v)
rk = RungeKutta(ls,ls,0.05,:BE_1_0_1)
opRK = TransientRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
sol_rk = solve(rk,opRK,u0,t0,T)
for (uh,t) in sol_rk
  println(uh)
end


createpvd("my_tests/transient_sol/rk_og") do pvd
  for (uh,t) in sol_rk
    pvd[t] = createvtk(Ω,"my_tests/transient_sol/rk_og_$t"*".vtu",cellfields=["u"=>uh])
  end
end
