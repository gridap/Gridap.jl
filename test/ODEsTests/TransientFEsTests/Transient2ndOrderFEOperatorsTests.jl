module Transient2nOrderFEOperatorsTests

using Gridap
using Test

# Analytical functions
u(x,t) = (1.0-x[1])*x[1]*(t^2+3.0)
u(t::Real) = x -> u(x,t)
v(t::Real) = ∂t(u)(t)
a(t::Real) = ∂tt(u)(t)
f(t) = x -> ∂tt(u)(x,t) + ∂t(u)(x,t) - Δ(u(t))(x)

u_const(x,t) = (1.0-x[1])*x[1]*(3.0)
u_const(t::Real) = x -> u_const(x,t)
v_const(t::Real) = ∂t(u_const)(t)
a_const(t::Real) = ∂tt(u_const)(t)
f_const(t) = x -> ∂tt(u_const)(x,t) + ∂t(u_const)(x,t) - Δ(u_const(t))(x)

domain = (0,1)
partition = (2,)
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

m(utt,v) = ∫(v*utt)dΩ
c(ut,v) = ∫(v*ut)dΩ
a(u,v) = ∫(∇(v)⊙∇(u))dΩ
b(t,v) = ∫(v*f(t))dΩ
b_const(v) = ∫(v*f_const(0.0))dΩ
m(t,utt,v) = m(utt,v)
c(t,ut,v) = c(ut,v)
a(t,u,v) = a(u,v)

res(t,u,v) = m(∂tt(u),v) + c(∂t(u),v) + a(u,v) - b(t,v)
jac(t,u,du,v) = a(du,v)
jac_t(t,u,dut,v) = c(dut,v)
jac_tt(t,u,dutt,v) = m(dutt,v)

op = TransientFEOperator(res,jac,jac_t,jac_tt,U,V0)
op_affine = TransientAffineFEOperator(m,c,a,b,U,V0)
op_const = TransientConstantFEOperator(m,c,a,b_const,U,V0)
op_const_mat = TransientConstantMatrixFEOperator(m,c,a,b,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1
γ = 0.5
β = 0.25

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)
vh0 = interpolate_everywhere(v(0.0),U0)
ah0 = interpolate_everywhere(a(0.0),U0)
vh0_const = interpolate_everywhere(v_const(0.0),U0)
ah0_const = interpolate_everywhere(a_const(0.0),U0)

ls = LUSolver()
ode_solver = Newmark(ls,dt,γ,β)

sol_t = solve(ode_solver,op,(uh0,vh0,ah0),t0,tF)
sol_affine_t = solve(ode_solver,op_affine,(uh0,vh0,ah0),t0,tF)
sol_const_t = solve(ode_solver,op_const,(uh0,vh0_const,ah0_const),t0,tF)
sol_const_mat_t = solve(ode_solver,op_const_mat,(uh0,vh0,ah0),t0,tF)

l2(w) = w*w

tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  @test tn≈_t_n
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

all_sol = [ (copy(uh_tn), tn) for (uh_tn, tn) in sol_t ]
all_el2 = [ sqrt(sum( ∫(l2( u(tn) - uhc_tn ))dΩ )) for (uhc_tn,tn) in all_sol ]
@test all( all_el2 .< tol )

_t_n = t0
for (uh_tn, tn) in sol_affine_t
  global _t_n
  _t_n += dt
  @test tn≈_t_n
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

all_sol = [ (copy(uh_tn), tn) for (uh_tn, tn) in sol_affine_t ]
all_el2 = [ sqrt(sum( ∫(l2( u(tn) - uhc_tn ))dΩ )) for (uhc_tn,tn) in all_sol ]
@test all( all_el2 .< tol )

_t_n = t0
for (uh_tn, tn) in sol_const_t
  global _t_n
  _t_n += dt
  @test tn≈_t_n
  e = u_const(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

all_sol = [ (copy(uh_tn), tn) for (uh_tn, tn) in sol_const_t ]
all_el2 = [ sqrt(sum( ∫(l2( u_const(tn) - uhc_tn ))dΩ )) for (uhc_tn,tn) in all_sol ]
@test all( all_el2 .< tol )

_t_n = t0
for (uh_tn, tn) in sol_const_mat_t
  global _t_n
  _t_n += dt
  @test tn≈_t_n
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
 @test el2 < tol
end

all_sol = [ (copy(uh_tn), tn) for (uh_tn, tn) in sol_const_mat_t ]
all_el2 = [ sqrt(sum( ∫(l2( u(tn) - uhc_tn ))dΩ )) for (uhc_tn,tn) in all_sol ]
@test all( all_el2 .< tol )

end
