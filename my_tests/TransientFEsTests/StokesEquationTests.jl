using Pkg
using Test
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap

u(x,t) = VectorValue(x[1],x[2])*t
u(t::Real) = x -> u(x,t)

p(x,t) = (x[1]-x[2])*t
p(t::Real) = x -> p(x,t)
q(x) = t -> p(x,t)

f(t) = x -> ∂t(u)(t)(x)-Δ(u(t))(x)+ ∇(p(t))(x)
g(t) = x -> (∇⋅u(t))(x)

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 2

reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
V0 = FESpace(
  model,
  reffeᵤ,
  conformity=:H1,
  dirichlet_tags="boundary"
)

reffeₚ = ReferenceFE(lagrangian,Float64,order-1)
Q = TestFESpace(
  model,
  reffeₚ,
  conformity=:H1,
  constraint=:zeromean
)

U = TransientTrialFESpace(V0,u)

P = TrialFESpace(Q)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

#
a(u,v) = ∫(∇(u)⊙∇(v))dΩ
b((v,q),t) = ∫(v⋅f(t))dΩ + ∫(q*g(t))dΩ
m(ut,v) = ∫(ut⋅v)dΩ

X = TransientMultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V0,Q])

lhs(t,(u,p),(v,q)) = ∫(u⋅v)
rhs(t,(u,p),(v,q)) =  ∫(v⋅f(t))dΩ + ∫(q*g(t))dΩ - ∫(∇(u)⊙∇(v))dΩ + ∫((∇⋅v)*p)dΩ - ∫(q*(∇⋅u))dΩ
jac(t,(u,p),(du,dp),(v,q)) = a(du,v) - ∫((∇⋅v)*dp)dΩ + ∫(q*(∇⋅du))dΩ
jac_t(t,(u,p),(dut,dpt),(v,q)) = m(dut,v)

b((v,q)) = b((v,q),0.0)


U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)
ph0 = interpolate_everywhere(p(0.0),P0)
xh0 = interpolate_everywhere([uh0,ph0],X0)

op = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,X,Y)
t0 = 0.0
tF = 1.0
dt = 0.1

ls = LUSolver()
ode_solver = EXRungeKutta(ls,dt,:EX_FE_1_0_1)

sol_t = solve(ode_solver,op,xh0,t0,tF)

l2(w) = w⋅w


tol = 1.0e-6
_t_n = t0

result = Base.iterate(sol_t)

for (xh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  uh_tn = xh_tn[1]
  ph_tn = xh_tn[2]
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  e = p(tn) - ph_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end
