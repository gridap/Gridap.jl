module DGFEOperatorsTests

using Gridap

u(x) = x[1] + x[2]
f(x) = 0.0
g(x) = u(x)

L = 1.0
limits = (0.0, L, 0.0, L)
ncellx = 4
model = CartesianDiscreteModel(domain=limits, partition=(ncellx,ncellx))

h = L / ncellx

γ = 10

order = 2
fespace = DLagrangianFESpace(Float64,model,order)

V = TestFESpace(fespace)
U = TrialFESpace(fespace)

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

btrian = BoundaryTriangulation(model,"boundary")
bquad = CellQuadrature(btrian,order=2)
nb = NormalVector(btrian)

strian = SkeletonTriangulation(model,"interior") # TODO
squad = CellQuadrature(strian,order=2)
ns = NormalVector(strian)


v = FEBasis(V)

vb = restrict(v,btrian)

qb = coordinates(bquad)
@show evaluate((∇(vb)*nb).cellbasis,qb)

uh = zero(U)

uhb = restrict(uh,btrian)

@show evaluate(nb,qb)
@show evaluate(∇(uhb),qb)

vs = restrict(v,strian)

qs = coordinates(squad)
@show evaluate((jump(vs*ns)).cellbasis1,qs)

kkk

a_Ω(v,u) = inner(∇(v), ∇(u))
b_Ω(v) = inner(v,f)
t_Ω = AffineFETerm(a_Ω,b_Ω,trian,quad)

a_∂Ω(v,u) = inner(v, ∇(u)*nb ) - inner(∇(v)*nb, u) + (γ/h) * inner(v,u)
b_∂Ω(v) = (γ/h) * inner(v,g) - inner(∇(v)*nb, g)
t_∂Ω = AffineFETerm(a_∂Ω,b_∂Ω,btrian,bquad)

a_Γ(v,u) = inner( jump(v*ns), mean(∇(u)) ) - inner( mean(∇(v)), jump(u*ns) ) + (γ/h) * inner( jump(v*ns), jump(u*ns))
t_Γ = LinearFETerm(a_Γ,strian,squad)

op = LinearFEOperator(V,U,t_Ω,t_∂Ω,t_Γ)

uh = solve(op)

e = u - uh

l2(u) = inner(u,u)
h1(u) = a_Ω(u,u) + l2(u)

el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8

end # module
