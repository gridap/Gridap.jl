module StokesEquationTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.FESpaces: get_algebraic_operator


θ = 0.5

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

res(t,(u,p),(v,q)) = a(u,v) + m(∂t(u),v) - ∫((∇⋅v)*p)dΩ + ∫(q*(∇⋅u))dΩ - b((v,q),t)
jac(t,(u,p),(du,dp),(v,q)) = a(du,v) - ∫((∇⋅v)*dp)dΩ + ∫(q*(∇⋅du))dΩ
jac_t(t,(u,p),(dut,dpt),(v,q)) = m(dut,v)

b((v,q)) = b((v,q),0.0)

mat((du1,du2),(v1,v2)) = a(du1,v1)+a(du2,v2)

U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)
ph0 = interpolate_everywhere(p(0.0),P0)
xh0 = interpolate_everywhere([uh0,ph0],X0)

op = TransientFEOperator(res,jac,jac_t,X,Y)

t0 = 0.0
tF = 1.0
dt = 0.1

ls = LUSolver()
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
  ph_tn = xh_tn[2]
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
  e = p(tn) - ph_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

all_sol  = [ (copy(xh_tn), tn) for (xh_tn, tn) in sol_t ]
all_el2u = [ sqrt(sum( ∫(l2( u(tn) - xhc_tn[1] ))dΩ )) for (xhc_tn,tn) in all_sol ]
all_el2p = [ sqrt(sum( ∫(l2( p(tn) - xhc_tn[2] ))dΩ )) for (xhc_tn,tn) in all_sol ]
@test all( all_el2u .< tol )
@test all( all_el2p .< tol )

end #module
