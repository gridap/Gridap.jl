module DGHeatEquationTests

using Gridap
using ForwardDiff
using LinearAlgebra
using Test
using Gridap.FESpaces: get_algebraic_operator

θ = 0.2

u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
u(t::Real) = x -> u(x,t)
f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)

L= 1.0
n = 2
domain = (0,L,0,L)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

order = 2

reffe = ReferenceFE(lagrangian,Float64,order)
V0 = FESpace(
  model,
  reffe,
  conformity=:L2
)
U = TransientTrialFESpace(V0)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,degree)
nb = get_normal_vector(Γ)

Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ,degree)
ns = get_normal_vector(Λ)

a(u,v) = ∫(∇(v)⋅∇(u))dΩ
b(v,t) = ∫(v*f(t))dΩ
m(u,v) = ∫(v*u)dΩ

h = 1.0 / n
γ = order*(order+1)
a_Γ(u,v) = ∫( (γ/h)*v*u - v*(∇(u)⋅nb) - (∇(v)⋅nb)*u )dΓ
b_Γ(v,t) = ∫( (γ/h)*v*u(t) - (∇(v)⋅nb)*u(t) )dΓ

a_Λ(u,v) = ∫( (γ/h)*jump(v*ns)⊙jump(u*ns) - jump(v*ns)⊙mean(∇(u)) - mean(∇(v))⊙jump(u*ns) )dΛ


res(t,u,v) = a(u,v) + m(∂t(u),v) + a_Γ(u,v) + a_Λ(u,v) - b(v,t) - b_Γ(v,t)
jac(t,u,du,v) = a(du,v) + a_Γ(du,v) + a_Λ(du,v)
jac_t(t,u,dut,v) = m(dut,v)

op = TransientFEOperator(res,jac,jac_t,U,V0)

t0 = 0.0
tF = 1.0
dt = 0.1

U0 = U(0.0)
uh0 = interpolate_everywhere(u(0.0),U0)

ls = LUSolver()
using Gridap.Algebra: NewtonRaphsonSolver
ode_solver = ThetaMethod(ls,dt,θ)

sol_t = solve(ode_solver,op,uh0,t0,tF)

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

all_sol = [ (copy(uh_tn), tn) for (uh_tn, tn) in sol_t ]
all_el2 = [ sqrt(sum( ∫(l2( u(tn) - uhc_tn ))dΩ )) for (uhc_tn,tn) in all_sol ]
@test all( all_el2 .< tol )

end #module
