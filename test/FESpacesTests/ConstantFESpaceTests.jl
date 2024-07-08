module ConstantFESpacesTests

using Gridap
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

Λ2=ConstantFESpace(model,field_type=VectorValue{2,Float64})
Gridap.FESpaces.test_fe_space(Λ2)
M2=TrialFESpace(Λ2)
a2(μ,λ) = ∫(λ⋅μ)dΩ
l2(λ) = ∫(VectorValue(0.0,0.0)⋅λ)dΩ
op2 = AffineFEOperator(a2,l2,M2,Λ2)
μ2h = solve(op2)
@assert sum(∫(μ2h⋅μ2h)dΩ) < 1.0e-12

end # module
