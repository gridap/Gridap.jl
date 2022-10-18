module TransientFEOperatorsTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.ODEs.ODETools
using Gridap.ODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator

θ = 0.4

# Analytical functions
u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*(t+3.0)
u(t::Real) = x -> u(x,t)
v(x) = t -> u(x,t)

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
  dirichlet_tags="boundary")

U = TransientTrialFESpace(V0,u)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

a(u,v) = ∫(∇(v)⊙∇(u))dΩ
m(u,v) = ∫(v*u)dΩ
b(v,t) = ∫(v*f(t))dΩ

res(t,u,v) = a(u,v) + m(∂t(u),v) - b(v,t)
jac(t,u,du,v) = a(du,v)
jac_t(t,u,dut,v) = m(dut,v)

op = TransientFEOperator(res,jac,jac_t,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1


U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ls = LUSolver()
ode_solver = ThetaMethod(ls,dt,θ)

sol_t = solve(ode_solver,op,uh0,t0,tF)

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

@test length( [uht for uht in sol_t] ) == ceil((tF - t0)/dt)

#

u0 = get_free_dof_values(uh0)
uf = get_free_dof_values(uh0)

odeop = get_algebraic_operator(op)

ode_cache = allocate_cache(odeop)
vθ = similar(u0)
nl_cache = nothing

# tf = t0+dt

ode_solver.θ == 0.0 ? dtθ = dt : dtθ = dt*ode_solver.θ
tθ = t0+dtθ
ode_cache = update_cache!(ode_cache,odeop,tθ)

using Gridap.ODEs.ODETools: ThetaMethodNonlinearOperator
nlop = ThetaMethodNonlinearOperator(odeop,tθ,dtθ,u0,ode_cache,vθ)

nl_cache = solve!(uf,ode_solver.nls,nlop,nl_cache)

K = nl_cache.A
h = nl_cache.b

# Steady version of the problem to extract the Laplacian and mass matrices
# tf = 0.1
tf = tθ
Utf = U(tf)
# fst(x) = -Δ(u(tf))(x)
fst(x) = f(tf)(x)
a(u,v) = ∫(∇(v)⊙∇(u))dΩ

function extract_matrix_vector(a,fst)
  btf(v) = ∫(v*fst)dΩ
  op = AffineFEOperator(a,btf,Utf,V0)
  ls = LUSolver()
  solver = LinearFESolver(ls)
  uh = solve(solver,op)

  tol = 1.0e-6
  e = uh-u(tf)
  l2(e) = e*e
  l2e = sqrt(sum( ∫(l2(e))dΩ ))
  # @test l2e < tol

  Ast = op.op.matrix
  bst = op.op.vector

  @test uh.free_values ≈ Ast \ bst

  return Ast, bst
end

A,rhs = extract_matrix_vector(a,fst)

gst(x) = u(tf)(x)
m(u,v) = ∫(u*v)dΩ

M,_ = extract_matrix_vector(m,gst)

@test rhs ≈ h
@test A+M/(θ*dt) ≈ K

rhs
h


end #module
