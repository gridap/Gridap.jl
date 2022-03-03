module StokesEquationAutoDiffTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.ODEs.ODETools
using Gridap.ODEs.TransientFETools
using Gridap.FESpaces
using Gridap.Arrays: test_array

# using Gridap.ODEs.ODETools: ThetaMethodLinear
import Gridap: ∇
import Gridap.ODEs.TransientFETools: ∂t

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

X₀ = evaluate(X,nothing)
dy = get_fe_basis(Y)
dx = get_trial_fe_basis(X₀)
xh = FEFunction(X₀,rand(num_free_dofs(X₀)))
xh_t = TransientCellField(xh,(xh,))

cell_j = get_array(jac(0.5,xh_t,dx,dy))
cell_j_t = get_array(jac_t(0.5,xh_t,dx,dy))

cell_j_auto = get_array(jacobian(x->res(0.5,TransientCellField(x,(xh,)),dy),xh))
cell_j_t_auto = get_array(jacobian(x->res(0.5,TransientCellField(xh,(x,)),dy),xh))

for i in 1:length(cell_j)
  test_array(cell_j[i].array[1,1],cell_j_auto[i].array[1,1],≈)
  test_array(cell_j[i].array[1,2],cell_j_auto[i].array[1,2],≈)
  test_array(cell_j[i].array[2,1],cell_j_auto[i].array[2,1],≈)
  test_array(cell_j_t[i].array[1,1],cell_j_t_auto[i].array[1,1],≈)
end

op = TransientFEOperator(res,X,Y)

U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)
ph0 = interpolate_everywhere(p(0.0),P0)
xh0 = interpolate_everywhere([uh0,ph0],X0)

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
  e = p(tn) - ph_tn
  el2 = sqrt(sum( ∫(l2(e))dΩ ))
  @test el2 < tol
end

end #module
