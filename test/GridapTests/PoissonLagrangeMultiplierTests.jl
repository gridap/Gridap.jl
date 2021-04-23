module PoissonLagrangeMultiplierTests

using Test
using Gridap
using Gridap.Geometry

function is_left(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  ( x[1] <= 0.5 )
end

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
bgface_to_mask = get_face_mask(labels,"boundary",1)
Γface_to_bgface = findall(bgface_to_mask)
model_Γ = BoundaryDiscreteModel(Polytope{1},model,Γface_to_bgface)

Γ = Triangulation(model_Γ)
Γface_coords = get_cell_coordinates(Γ)
Γface_mask = lazy_map(is_left,Γface_coords)
Γlface_Γface = findall(Γface_mask)
Γrface_Γface = findall(!,Γface_mask)
Γl = BoundaryTriangulation(model,view(Γface_to_bgface,Γlface_Γface))
Γr = BoundaryTriangulation(model,view(Γface_to_bgface,Γrface_Γface))

Λ = SkeletonTriangulation(model_Γ)

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
dΩ = Measure(Ω,degree)
dΓl = Measure(Γl,degree)
dΓr = Measure(Γr,degree)
dΛ = Measure(Λ,degree)
n = get_normal_vector(Γ)
nΛ = get_normal_vector(Λ)

uₑ(x) = x[1]^2 + x[2]^2
f(x) = -Δ(uₑ)(x)

# Weak form. Additional non needed terms are added for testing purposes
a((u,λ),(v,η)) = ∫( ∇(u)⋅∇(v) )dΩ +
               ∫( (u+λ)*(η+v) - u*v - η*λ + 0.0*(n⋅∇(u)-λ)*(n⋅∇(v)-η) )dΓl +
               ∫( (u+λ)*(η+v) - u*v - η*λ + 0.0*(n⋅∇(u)-λ)*(n⋅∇(v)-η) )dΓr +
               ∫( 0.0*jump(n⋅∇(u))*mean(η) )dΛ +
               ∫( 0.0*jump(nΛ⋅∇(u))*mean(η) )dΛ
l((v,η)) = ∫( f*v )dΩ + ∫( uₑ*η )dΓl + ∫( uₑ*η )dΓr

op = AffineFEOperator(a,l,X,Y)
uh,λh = solve(op)

e = uh-uₑ
l2(v) = √(∑(∫(v*v)dΩ))

tol = 1.0e-10
@test l2(e) < tol

end
