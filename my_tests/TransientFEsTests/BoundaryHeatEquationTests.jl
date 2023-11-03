using Pkg
using Test
using ForwardDiff
using LinearAlgebra
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap


# Analytical functions
# u(x,t) = (x[1]+x[2])*t
# u(x,t) = (2*x[1]+x[2])*t
# u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(x,t) = (1.0-x[1])*x[1]*t

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
b_Γ(v,t) = ∫(v*(∇(u(t))⋅nb))dΓ

lhs(t,u,v) = ∫( v* (u) )dΩ
rhs(t,u,v) = b(v,t) + b_Γ(v,t) - a(u,v)

lhs(t,u,v) = ∫( v* (u) )dΩ
rhs(t,u,v) = ∫(v*f(t))dΩ -  ∫(( ∇(v)⊙∇(u) ))dΩ
jac(t,u,du,v) = ∫(( ∇(v)⊙∇(du) ))dΩ
jac_t(t,u,dut,v) = ∫( dut*v )dΩ

op = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ls = LUSolver()

ode_solver = EXRungeKutta(ls,dt,:EX_FE_1_0_1)

sol_t = solve(ode_solver,op,uh0,t0,tF)

l2(w) = w*w

tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  println(el2)
  # println(@test el2 < tol)
end

a(u,v) = ∫(∇(v)⋅∇(u))dΩ
b(v,t) = ∫(v*f(t))dΩ
m(ut,v) = ∫(ut*v)dΩ
b_Γ(v,t) = ∫(v*(∇(u(t))⋅nb))dΓ
res(t,u,v) = a(u,v) + m(∂t(u),v) - b(v,t) - b_Γ(v,t)
op = TransientFEOperator(res,jac,jac_t,U,V0)
ode_solver = ThetaMethod(ls,dt,0.2)
