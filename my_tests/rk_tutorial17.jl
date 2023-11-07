using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap
# using LaTeXStrings
using Plots
# using JLD


n = 16
p = 2
degree = 4*(p+1)
L = 1
dx = L/n
dt = 0.001
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

u0 = interpolate_everywhere(0.0,U(0.0))


ls = LUSolver()

m(t,u,v) = ∫(v*u)dΩ

f(t) = sin(π*t)
lhs(t,u,v) = ∫( v* (u) )dΩ
rhs(t,u,v) = ∫(v*f(t))dΩ -  ∫(( ∇(v)⊙∇(u) ))dΩ
jac(t,u,du,v) = ∫((  ∇(v)⊙∇(du) ))dΩ
jac_t(t,u,dut,v) = ∫( dut*v )dΩ


#### Solve with standard EXRungeKutta with FE table
rk_fe = EXRungeKutta(ls,dt,:EX_FE_1_0_1)
# rk_fe = EXRungeKutta(ls,0.0001,:EX_SSP_3_0_3)
opRK_fe = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V)
sol_rk_fe = solve(rk_fe,opRK_fe,u0,t0,T)


createpvd("my_tests/transient_sol/poisson_transient_solution") do pvd
  for (uh,t) in sol_rk_fe
    pvd[t] = createvtk(Ω,"my_tests/transient_sol/poisson_transient_solution_$t"*".vtu",cellfields=["u"=>uh])
  end
end
