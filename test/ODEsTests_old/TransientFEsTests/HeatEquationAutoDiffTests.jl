module HeatEquationAutoDiffTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.ODEs.ODETools
using Gridap.ODEs.TransientFETools
using Gridap.FESpaces: get_algebraic_operator
using Gridap.Arrays: test_array

θ = 0.2

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t::Real) = x -> u(x,t)
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
  dirichlet_tags="boundary"
)
U = TransientTrialFESpace(V0,u)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

#
a(u,v) = ∫(∇(v)⋅∇(u))dΩ
b(v,t) = ∫(v*f(t))dΩ

res(t,u,v) = a(u,v) + ∫(∂t(u)*v)dΩ - b(v,t)
jac(t,u,du,v) = a(du,v)
jac_t(t,u,dut,v) = ∫(dut*v)dΩ

U₀ = evaluate(U,nothing)
dv = get_fe_basis(V0)
du = get_trial_fe_basis(U₀)
uh = FEFunction(U₀,rand(num_free_dofs(U₀)))
uh_t = TransientCellField(uh,(uh,))

cell_j = get_array(jac(0.5,uh_t,du,dv))
cell_j_t = get_array(jac_t(0.5,uh_t,du,dv))

cell_j_auto = get_array(jacobian(x->res(0.5,TransientCellField(x,(uh,)),dv),uh))
cell_j_t_auto = get_array(jacobian(x->res(0.5,TransientCellField(uh,(x,)),dv),uh))

test_array(cell_j_auto,cell_j,≈)
test_array(cell_j_t_auto,cell_j_t,≈)

op = TransientFEOperator(res,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ls = LUSolver()
using Gridap.Algebra: NewtonRaphsonSolver
nls = NLSolver(ls;show_trace=false,method=:newton) #linesearch=BackTracking())
ode_solver = ThetaMethod(ls,dt,θ)

sol_t = solve(ode_solver,op,uh0,t0,tF)

# Juno.@enter Base.iterate(sol_t)

l2(w) = w*w

tol = 1.0e-6
_t_n = t0

for (uh_tn, tn) in sol_t
  global _t_n
  _t_n += dt
  e = u(tn) - uh_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

end #module
