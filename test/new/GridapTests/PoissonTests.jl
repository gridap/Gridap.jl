module PoissonTests

using Test
using Gridap
import Gridap: ∇

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)
order = 2
T = Float64
diritags = "boundary"

u(x) = x[1]^2 + x[2]
∇u(x) = VectorValue( 2*x[1], one(x[2]) )
Δu(x) = 2
f(x) = - Δu(x)

∇(::typeof(u)) = ∇u

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

V = FESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 valuetype=T,
 dirichlet_tags=diritags)

U = TrialFESpace(V,u)

a(v,u) = inner(∇(v),∇(u))
l(v) = v*f
t_Ω = AffineFETerm(a,l,trian,quad)

op = AffineFEOperator(V,U,t_Ω)

uh = solve(op)

e = u - uh

l2(u) = inner(u,u)
sh1(u) = a(u,u)
h1(u) = sh1(u) + l2(u)

el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
ul2 = sqrt(sum( integrate(l2(uh),trian,quad) ))
uh1 = sqrt(sum( integrate(h1(uh),trian,quad) ))

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
