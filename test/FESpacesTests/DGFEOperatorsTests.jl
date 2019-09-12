module DGFEOperatorsTests

using Test
using Gridap
import Gridap: ∇

u(x) = x[1] + x[2]
∇u(x) = VectorValue(1.0,1.0)
∇(::typeof(u)) = ∇u
f(x) = 0.0
g(x) = u(x)

L = 1.0
limits = (0.0, L, 0.0, L)
ncellx = 4
model = CartesianDiscreteModel(domain=limits, partition=(ncellx,ncellx))

h = L / ncellx

γ = 10

order = 3
fespace = DLagrangianFESpace(Float64,model,order)

V = TestFESpace(fespace)
U = TrialFESpace(fespace)

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2*order)

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,order=2*order)
nb = NormalVector(btrian)

strian = SkeletonTriangulation(model)
squad = CellQuadrature(strian,order=2*order)
ns = NormalVector(strian)

a_Ω(v,u) = inner(∇(v), ∇(u))
b_Ω(v) = inner(v,f)
t_Ω = AffineFETerm(a_Ω,b_Ω,trian,quad)

a_∂Ω(v,u) = (γ/h) * inner(v,u) - inner(v, ∇(u)*nb ) - inner(∇(v)*nb, u)
b_∂Ω(v) = (γ/h) * inner(v,g) - inner(∇(v)*nb, g)
t_∂Ω = AffineFETerm(a_∂Ω,b_∂Ω,btrian,bquad)

a_Γ(v,u) = (γ/h) * inner( jump(v*ns), jump(u*ns)) - inner( jump(v*ns), mean(∇(u)) ) - inner( mean(∇(v)), jump(u*ns) ) 
t_Γ = LinearFETerm(a_Γ,strian,squad)

op = LinearFEOperator(V,U,t_Ω,t_∂Ω,t_Γ)

uh = solve(op)

e = u - uh

#writevtk(trian,"trian",cellfields=["uh"=>uh,"e"=>e])

l2(u) = inner(u,u)
h1(u) = a_Ω(u,u) + l2(u)

el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8

end # module
