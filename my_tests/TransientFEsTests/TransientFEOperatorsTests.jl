using Pkg
using Test
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap

# Analytical functions
u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*(t)
u(t::Real) = x -> u(x,t)
v(x) = t -> u(x,t)
f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)
∂tu(x,t) = ∂t(u)(x,t)
∂tu(t::Real) = x -> ∂tu(x,t)

# Domain and triangulations
domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = FESpace(
  model,
  reffe,
  conformity=:H1,
  dirichlet_tags="boundary")
U = TransientTrialFESpace(V0,u)
Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

# Affine FE operator
lhs(t,u,v) = ∫( v* (u) )dΩ
rhs(t,u,v) = ∫(v*f(t))dΩ -  ∫(( ∇(v)⊙∇(u) ))dΩ
jac(t,u,du,v) = ∫(( ∇(v)⊙∇(du) ))dΩ
jac_t(t,u,dut,v) = ∫( dut*v )dΩ


op = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V0)

ls = LUSolver()

# Time stepping
t0 = 0.0
tF = 1.0
dt = 0.01

# Initial solution
U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ode_solver = EXRungeKutta(ls,dt,:EX_FE_1_0_1)

sol_t = solve(ode_solver,op,uh0,t0,tF)

l2(w) = w*w

tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  println(el2)
  @test el2 < tol
end
