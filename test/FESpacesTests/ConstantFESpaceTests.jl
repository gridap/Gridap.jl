
module ConstantFESpacesTests

using Gridap
using GridapBiotElasticity
using Test

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)
Λ=ConstantFESpace(model)
Gridap.FESpaces.test_fe_space(Λ)
M=TrialFESpace(Λ)

order = 2
u((x,y)) = (x+y)^order
f(x) = -Δ(u,x)
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(u,V)
Y = MultiFieldFESpace([V,Λ])
X = MultiFieldFESpace([U,M])

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
a((u,μ),(v,λ)) = ∫(∇(v)⋅∇(u))dΩ + ∫(v*μ)dΩ + ∫(λ*u)dΩ
l((v,λ)) = ∫(λ*u)dΩ + ∫(v*f)dΩ
op = AffineFEOperator(a,l,X,Y)
uh = solve(op)

@assert sum(∫((uh[1]-u)*(uh[1]-u))dΩ) < 1.0e-14
abs(sum(∫(uh[2])dΩ)) < 1.0e-12

end # module
