module TransientFETests

using Gridap
using Test
using Gridap.ODEs.ODETools
using Gridap.ODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

u(x,t) = (x[1] + x[2])*t
u(t::Real) = x -> u(x,t)
∇u(x,t) = VectorValue(1,1)*t
∇u(t::Real) = x -> ∇u(x,t)
import Gridap: ∇
∇(::typeof(u)) = ∇u
∇(u) === ∇u

θ = 1.0

∂tu(t) = x -> x[1]+x[2]
import Gridap.ODEs.TransientFETools: ∂t
∂t(::typeof(u)) = ∂tu
@test ∂t(u) === ∂tu

f(t) = x -> (x[1]+x[2])

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = TestFESpace(
  model,
  reffe,
  conformity=:H1,
  dirichlet_tags="boundary"
)

U = TransientTrialFESpace(V0,u)
U0 = TrialFESpace(V0,u(0.0))
@test test_transient_trial_fe_space(U)

U0 = U(1.0)
ud0 = copy(get_dirichlet_dof_values(U0))
_ud0 = get_dirichlet_dof_values(U0)
U1 = U(2.0)
ud1 = copy(get_dirichlet_dof_values(U1))
_ud1 = get_dirichlet_dof_values(U1)
@test all(ud0 .≈ 0.5ud1)
all(_ud0 .≈ _ud1)

Ut = ∂t(U)
Ut.dirichlet_t
Ut0 = Ut(0.0)
Ut0.dirichlet_values

Ut1 = Ut(1.0)
utd0 = copy(get_dirichlet_dof_values(Ut0))
utd1 = copy(get_dirichlet_dof_values(Ut1))
@test all(utd0 .== utd1)
@test all(utd1 .== ud0)

Ω = Triangulation(model)
degree = 2
dΩ = Measure(Ω,degree)

a(u,v) = ∫(∇(v)⋅∇(u))dΩ
b(v,t) = ∫(v*f(t))dΩ

res(t,u,v) = a(u,v) + ∫(∂t(u)*v)dΩ - b(v,t)
jac(t,u,du,v) = a(du,v)
jac_t(t,u,dut,v) = ∫(dut*v)dΩ

U0 = U(0.0)
_res(u,v) = a(u,v) + 10.0*∫(u*v)dΩ - b(v,0.0)
_jac(u,du,v) = a(du,v) + 10.0*∫(du*v)dΩ
_op = FEOperator(_res,_jac,U0,V0)

uh = interpolate_everywhere(0.0,U0)#1.0)
using Gridap.FESpaces: allocate_residual, allocate_jacobian
_r = allocate_residual(_op,uh)
_J = allocate_jacobian(_op,uh)
using Gridap.FESpaces: residual!, jacobian!
residual!(_r,_op,uh)
jacobian!(_J,_op,uh)

op = TransientFEOperator(res,jac,jac_t,U,V0)
odeop = get_algebraic_operator(op)
cache = allocate_cache(odeop)

r = allocate_residual(op,uh,cache)
J = allocate_jacobian(op,uh,cache)
uh10 = interpolate_everywhere(0.0,U0)#10.0)
xh = TransientCellField(uh,(uh10,))
residual!(r,op,0.0,xh,cache)
jacobian!(J,op,1.0,xh,1,1.0,cache)
jacobian!(J,op,1.0,xh,2,10.0,cache)
@test all(r.≈_r)
@test all(J.≈_J)

U0 = U(0.0)
uh0 = interpolate_everywhere(0.0,U0)
@test test_transient_fe_operator(op,uh0)

u0 = u(0.0)
t0 = 0.0
tF = 1.0
dt = 0.1

ls = LUSolver()
# using LineSearches: BackTracking
tol = 1.0
maxiters = 20
using Gridap.Algebra: NewtonRaphsonSolver
nls = NLSolver(ls;show_trace=true,method=:newton) #linesearch=BackTracking())
ode_solver = ThetaMethod(nls,dt,1.0)
@test test_transient_fe_solver(ode_solver,op,uh0,t0,tF)

xh = TransientCellField(uh,(uh,))
residual!(r,op,0.1,xh,cache)
jacobian!(J,op,1.0,xh,1,1.0,cache)
jacobian!(J,op,1.0,xh,2,10.0,cache)

u0 = get_free_dof_values(uh0)
solver = ode_solver
t0 = 0.0
ode_cache = allocate_cache(odeop)
cache = nothing
uf = copy(u0)
dt = solver.dt
tf = t0+dt
update_cache!(ode_cache,odeop,tf)
using Gridap.ODEs.ODETools: ThetaMethodNonlinearOperator
vf = copy(u0)
nlop = ThetaMethodNonlinearOperator(odeop,tf,dt,u0,ode_cache,vf)

x = copy(nlop.u0)

b1 = allocate_residual(nlop,x)
residual!(b1,nlop,x)
b2 = allocate_residual(nlop,x)
residual!(b2,nlop.odeop,nlop.tθ,(x,10.0*x),nlop.ode_cache)
@test all(b1 .≈ b2)
J1 = allocate_jacobian(nlop,x)
jacobian!(J1,nlop,x)
J2 = allocate_jacobian(nlop,x)
jacobian!(J2,nlop.odeop,nlop.tθ,(x,10.0*x),1,1.0,nlop.ode_cache)
jacobian!(J2,nlop.odeop,nlop.tθ,(x,10.0*x),2,10.0,nlop.ode_cache)
@test all(J1 .≈ J2)
using Gridap.Algebra: test_nonlinear_operator
test_nonlinear_operator(nlop,x,b1,jac=J1)

x .= 0.0
r = allocate_residual(nlop,x)
residual!(r,nlop,x)
J = allocate_jacobian(nlop,x)
jacobian!(J,nlop,x)

cache = solve!(uf,solver.nls,nlop)
df = cache.df
ns = cache.ns

function linsolve!(x,A,b)
  numerical_setup!(ns,A)
  solve!(x,ns,b)
end

p = copy(x)
p .= 0.0
l_sol = linsolve!(p,J,-r)
J*l_sol .≈ -r
x = x + l_sol
@test all(abs.(residual!(r,nlop,x)) .< 1e-6)

residual!(r,nlop,x)
jacobian!(J,nlop,x)
p .= 0.0
l_sol = linsolve!(p,J,-r)

cache = solve!(uf,solver.nls,nlop)
@test all(uf .≈ x)
solve!(uf,solver.nls,nlop,cache)
@test all(uf .≈ x)

uf .= 0.0
x = copy(nlop.u0)
cache = Gridap.Algebra._new_nlsolve_cache(x,nls,nlop)
df = cache.df
ns = cache.ns
x .= 0.0
l_sol = linsolve!(x,df.DF,df.F)
@test all(df.DF*l_sol.≈df.F)
x .= 0
Gridap.Algebra.nlsolve(df,x;linsolve=linsolve!,nls.kwargs...)

using Gridap.FESpaces: get_algebraic_operator
odeop = get_algebraic_operator(op)
sol_ode_t = solve(ode_solver,odeop,u0,t0,tF)

test_ode_solution(sol_ode_t)
_t_n = t0
for (u_n, t_n) in sol_ode_t
  global _t_n
  _t_n += dt
  @test t_n≈_t_n
  @test all(u_n .≈ t_n)
end

ode_solver = ThetaMethod(nls,dt,θ)
sol_ode_t = solve(ode_solver,odeop,u0,t0,tF)
test_ode_solution(sol_ode_t)
_t_n = t0
un, tn = Base.iterate(sol_ode_t)
for (u_n, t_n) in sol_ode_t
  global _t_n
  _t_n += dt
  @test t_n≈_t_n
  @test all(u_n .≈ t_n)
end

sol_t = solve(ode_solver,op,uh0,t0,tF)
@test test_transient_fe_solution(sol_t)

_t_n = 0.0
for (u_n, t_n) in sol_t
  global _t_n
  _t_n += dt
  @test t_n≈_t_n
  @test all(u_n.free_values .≈ t_n)
end

l2(w) = w*w

# h1(w) = a(w,w) + l2(w)

_t_n = t0
for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  @test tn≈_t_n
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
  # writevtk(trian,"sol at time: $tn",cellfields=["u" => uh_tn])
end

end #module
