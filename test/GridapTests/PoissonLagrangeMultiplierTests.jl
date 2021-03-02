module PoissonLagrangeMultiplierTests

using Test
using Gridap
using Gridap.Geometry

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
bgface_to_mask = get_face_mask(labels,"boundary",1)
model_Γ = BoundaryDiscreteModel(Polytope{1},model,bgface_to_mask)

order = 2
reffe_u = ReferenceFE(lagrangian,Float64,order)
reffe_λ = ReferenceFE(lagrangian,Float64,order-1)
V = TestFESpace(model,reffe_u,conformity=:H1)
S = TestFESpace(model_Γ,reffe_λ,conformity=:L2)
U = TrialFESpace(V)
L = TrialFESpace(S)
Y = MultiFieldFESpace([V,S])
X = MultiFieldFESpace([U,L])

degree = 2*order
Ω = Triangulation(model)
Γ = Triangulation(model_Γ)
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

uₑ(x) = x[1]^2 + x[2]^2
f(x) = -Δ(uₑ)(x)

a((u,λ),(v,η)) = ∫( ∇(u)⋅∇(v) )dΩ + ∫( u*η + λ*v )dΓ
l((v,η)) = ∫( f*v )dΩ + ∫( uₑ*η )dΓ

op = AffineFEOperator(a,l,X,Y)
uh,λh = solve(op)

e = uh-uₑ
l2(v) = √(∑(∫(v*v)dΩ))

tol = 1.0e-10
@test l2(e) < tol

end
