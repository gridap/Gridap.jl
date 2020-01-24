module PoissonDGTests

using Test
using Gridap
import Gridap: ∇

#domain = (0,1,0,1)
#partition = (4,4)
#model = CartesianDiscreteModel(domain,partition)
#const h = (domain[2]-domain[1]) / partition[1]

using Gridap.Geometry: DiscreteModelMock
model = DiscreteModelMock()
const h = 1

order = 2
const γ = 10


trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model)
bdegree = order
bquad = CellQuadrature(btrian,bdegree)
const bn = get_normal_vector(btrian)

strian = SkeletonTriangulation(model)
sdegree = order
squad = CellQuadrature(strian,sdegree)
const sn = get_normal_vector(strian)

u_scal(x) = x[1]^2 + x[2]
∇u_scal(x) = VectorValue( 2*x[1], one(x[2]) )
Δu_scal(x) = 2
f_scal(x) = - Δu_scal(x)
∇(::typeof(u_scal)) = ∇u_scal

T = Float64
u = u_scal
f = f_scal

V = TestFESpace(
 model=model,
 order=order,
 reffe=:Lagrangian,
 conformity=:L2,
 valuetype=T)

U = TrialFESpace(V,u)

a(v,u) = inner(∇(v),∇(u))
l(v) = v*f
t_Ω = AffineFETerm(a,l,trian,quad)

a_Γd(v,u) = (γ/h)*v*u  - v*(bn*∇(u)) - (bn*∇(v))*u
l_Γd(v) = (γ/h)*v*u - (bn*∇(v))*u
t_Γd = AffineFETerm(a_Γd,l_Γd,btrian,bquad)

a_Γ(v,u) = (γ/h)*jump(v*sn)*jump(u*sn) - jump(v*sn)*mean(∇(u)) -  mean(∇(v))*jump(u*sn)
t_Γ = LinearFETerm(a_Γ,strian,squad)

op = AffineFEOperator(V,U,t_Ω,t_Γ,t_Γd)

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

