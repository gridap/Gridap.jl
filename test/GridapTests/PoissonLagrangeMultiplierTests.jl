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

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(Γ)

D = num_cell_dims(Ω)
glue = get_glue(Γ,Val(D-1))
Γface_to_face = glue.tface_to_mface
Γface_to_coords = get_cell_coordinates(Γ)
Γface_mask = lazy_map(is_left,Γface_to_coords)
Γlface_Γface = findall(Γface_mask)
Γrface_Γface = findall(!,Γface_mask)
Γl = BoundaryTriangulation(model,view(Γface_to_face,Γlface_Γface))
Γr = BoundaryTriangulation(model,view(Γface_to_face,Γrface_Γface))

order = 2
reffe_u = ReferenceFE(lagrangian,Float64,order)
reffe_λ = ReferenceFE(lagrangian,Float64,order-1)
V = TestFESpace(Ω,reffe_u,conformity=:H1)
S = TestFESpace(Γ,reffe_λ,conformity=:L2)
U = TrialFESpace(V)
L = TrialFESpace(S)
Y = MultiFieldFESpace([V,S])
X = MultiFieldFESpace([U,L])

#uh, λh = FEFunction(Y,rand(num_free_dofs(Y)))
#writevtk(Ω,"t_Ω",nsubcells=10,cellfields=["uh"=>uh])
#writevtk(Γ,"t_Γ",nsubcells=10,cellfields=["uh"=>uh,"λh"=>λh])
#writevtk(Λ,"t_Λ",cellfields=["uh"=>mean(uh),"λh"=>mean(λh)])
#writevtk(Γr,"t_Γr",nsubcells=10,cellfields=["uh"=>uh,"λh"=>λh])
#writevtk(Γl,"t_Γl",nsubcells=10,cellfields=["uh"=>uh,"λh"=>λh])

degree = 2*order
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
