module PLaplacianWithAutodiffTests

using Test
using Gridap

u(x) = x[1] + x[2]
f(x) = 0
const p = 3
@law j(∇u) = norm(∇u)^(p-2) * ∇u

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 3
T = Float64

V = TestFESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 conformity=:H1,
 valuetype=T,
 dirichlet_tags="boundary")

U = TrialFESpace(V,u)

degree = 2*order
trian = get_triangulation(model)
quad = CellQuadrature(trian,degree)

res(u,v) = ∇(v)⋅j(∇(u)) - v*f

t_Ω = FETerm(res,trian,quad)

op = FEOperator(U,V,t_Ω)

nls = NLSolver(show_trace=false, method=:newton)
solver = FESolver(nls)

x = rand(T,num_free_dofs(U))
uh0 = FEFunction(U,x)
uh, = solve!(uh0,solver,op)

e = u - uh

el2 = sqrt(sum(integrate(e*e,quad)))
eh1 = sqrt(sum(integrate(e*e+∇(e)⋅∇(e),quad)))

@test el2 < 1.e-8
@test eh1 < 1.e-7

end # module
