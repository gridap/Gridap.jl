module EmptyDomainsTests

using Gridap
import Gridap: ∇
using Gridap.Geometry, Gridap.CellData
using Test

# Analytical functions
u(x) = VectorValue( x[1]^2 + 2*x[2]^2, -x[1]^2 )
∇u(x) = TensorValue( 2*x[1], 4*x[2], -2*x[1], zero(x[1]) )
Δu(x) = VectorValue( 6, -2 )
p(x) = x[1] + 3*x[2]
∇p(x) = VectorValue(1,3)
s(x) = -Δu(x)
f(x) = -Δu(x) + ∇p(x)
g(x) = tr(∇u(x))
∇(::typeof(u)) = ∇u
∇(::typeof(p)) = ∇p

# Mesh
model = CartesianDiscreteModel((0,1,0,1),(10,10))
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[6,7,8])

# Domains and measures
Ω = Interior(model,Int[])
Ωs = Interior(model,Int[])
Ωf = Interior(model,Int[])
Λ = Boundary(model,Int[])
Γ = Interface(Ωf,Ωs)
n_Λ = get_normal_vector(Λ)
n_Γ = get_normal_vector(Γ)
k = 2
degree = 2*k
dΩ = Measure(Ω,degree)
dΩs = Measure(Ωs,degree)
dΩf = Measure(Ωf,degree)
dΛ = Measure(Λ,degree)
dΓ = Measure(Γ,degree)

# FE Spaces
reffe_u = ReferenceFE(lagrangian,VectorValue{2,Float64},k)
reffe_p = ReferenceFE(lagrangian,Float64,k-1,space=:P)
V = TestFESpace(Ω,reffe_u,dirichlet_tags="dirichlet")
Q = TestFESpace(Ωf,reffe_p)
U = TrialFESpace(V,u)
P = Q
Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

# Weak form
a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) )*dΩs +
  ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )*dΩf

l((v,q)) =
  ∫( v⋅s )*dΩs +
  ∫( v⋅f + q*g )*dΩf +
  ∫( v⋅(n_Λ⋅∇u) - (n_Λ⋅v)*p )*dΛ +
  ∫( - (n_Γ.⁺⋅v.⁺)*p )*dΓ

# FE problem
op = AffineFEOperator(a,l,X,Y)
A = get_matrix(op)
b = get_vector(op)
@test size(A) == (0,0)
@test size(b) == (0,)
uh, ph = FEFunction(X,Float64[])

# Errors
eu = u - uh
ep = p - ph
eu_l2 = sqrt(sum(∫( eu⋅eu )*dΩ))
eu_h1 = sqrt(sum(∫( eu⋅eu + ∇(eu)⊙∇(eu) )*dΩ))
ep_l2 = sqrt(sum(∫( ep*ep )*dΩf))
tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

###############################################
# Multifield jumps

V = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1))
W = MultiFieldFESpace([V,V])

Λ = Skeleton(model,falses(num_cells(model)))
dΛ = Measure(Λ,2)

a((u,p),(v,q)) = ∫(jump(∇(u)) ⊙ jump(∇(v)))dΛ
assemble_matrix(a,W,W)

###############################################
# Multidomain

Ωa = Triangulation(model,collect(Int32,1:50))
Ωb = Triangulation(model,collect(Int32,51:100))
Ωc = Triangulation(model,Int32[])
Γ = Interface(Ωa,Ωb)

Va = TestFESpace(Ωa, ReferenceFE(lagrangian, Float64, 1))
Vb = TestFESpace(Ωb, ReferenceFE(lagrangian, Float64, 1))
X = MultiFieldFESpace([Va,Vb])

Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
dΩa = Measure(Ωa,2)
dΩb = Measure(Ωb,2)
dΩc = Measure(Ωc,2)
dΓ = Measure(Γ,2)
a1((ua,ub),(va,vb)) = ∫( ∇(va)⊙∇(ub) + vb*ub )*dΩa + ∫( ∇(vb)⊙∇(ua) + va*(ua+ub))*dΩb + ∫( mean(va*ub) + mean(vb*ua) )*dΓ
assemble_matrix(a1,X,X)
a2((ua,ub),(va,vb)) = ∫(Π(ua,Ωb)*Π(vb,Ωa))dΩc
assemble_matrix(a2,X,X)

###############################################
# Autodiff when the target domain is empty

V = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1))
uh = zero(V)
vh = get_fe_basis(V)

# Case 1
Ω_empty = Triangulation(model, Int32[])
dΩ = Measure(Ω_empty, 2)
j1(u) = ∫(u * u)dΩ
dj1 = gradient(j1, uh)
dj1_vec = assemble_vector(dj1,V)
j2(u) = ∫(u*vh)dΩ
dj2 = jacobian(j2, uh)
dj2_mat = assemble_matrix(dj2,V,V)

# Case 2
Ωa = Triangulation(model,Int32[])
Ωb = Triangulation(model,Int32[1])
Ω = AppendedTriangulation(Ωa,Ωb)
dΩ = Measure(Ω,2)
j3(u) = ∫(u)dΩ
dj3 = gradient(j3,uh)
dj3_vec = assemble_vector(dj3,V)

j4(u) = ∫(u*vh)dΩ
dj4 = jacobian(j4,uh)
dj4_mat = assemble_matrix(dj4,V,V)

end # module
