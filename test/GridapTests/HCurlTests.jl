module HCurlTests

using Gridap
using LinearAlgebra
using Test

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition) |> simplexify
order = 0
u((x,y)) = 2*VectorValue(-y,x)
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(V,u)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*(order+1))
a(u,v) = ∫( (∇×u)⋅(∇×v) + u⋅v )dΩ
l(v) = ∫(u⋅v)dΩ
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition) |> simplexify
order = 0
u((x,y,z)) = 2*VectorValue(-y,x,0.) - VectorValue(0.,-z,y)
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(V,u)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*(order+1))
a(u,v) = ∫( (∇×u)⋅(∇×v) + u⋅v )dΩ
l(v) = ∫(u⋅v)dΩ
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)
order = 1
u((x,y,z)) = 2*VectorValue(-y,x,0.) - VectorValue(0.,-z,y)
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(V,u)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*(order+1))
a(u,v) = ∫( (∇×u)⋅(∇×v) + u⋅v )dΩ
l(v) = ∫(u⋅v)dΩ
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

end # module
