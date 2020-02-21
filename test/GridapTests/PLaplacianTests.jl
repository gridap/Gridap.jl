module PLaplacianTests

using Test
using Gridap
import Gridap: ∇
using LinearAlgebra: norm
using Gridap.Geometry: DiscreteModelMock

model = DiscreteModelMock()

u(x) = x[1] + x[2]
∇u(x) = VectorValue(1,1)
f(x) = 0
∇(::typeof(u)) = ∇u

const p = 3
@law flux(∇u) = norm(∇u)^(p-2) * ∇u
@law dflux(∇du,∇u) = (p-2)*norm(∇u)^(p-4)*inner(∇u,∇du)*∇u + norm(∇u)^(p-2)*∇du

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

degree = order
trian = get_triangulation(model)
quad = CellQuadrature(trian,degree)

res(u,v) = inner( ∇(v), flux(∇(u)) ) - inner(v,f)
jac(u,v,du) = inner(  ∇(v) , dflux(∇(du),∇(u)) )

t_Ω = FETerm(res,jac,trian,quad)

op = FEOperator(V,U,t_Ω)

nls = NLSolver(show_trace=false, method=:newton)
solver = FESolver(nls)

x = rand(T,num_free_dofs(U))
uh0 = FEFunction(U,x)
uh, = solve!(uh0,solver,op)

e = u - uh

l2(u) = inner(u,u)
sh1(u) = inner(∇(u),∇(u))
h1(u) = sh1(u) + l2(u)

el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
ul2 = sqrt(sum( integrate(l2(uh),trian,quad) ))
uh1 = sqrt(sum( integrate(h1(uh),trian,quad) ))

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7


end # module
