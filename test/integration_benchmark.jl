
using BenchmarkTools, ProfileView
using SparseArrays, LinearAlgebra
using Gridap, Gridap.FESpaces, Gridap.Algebra

model = CartesianDiscreteModel((0,1,0,1),(100,100))

reffe_u = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
reffe_p = ReferenceFE(lagrangian, Float64, 1)
reffe_j = ReferenceFE(raviart_thomas, Float64, 0)
U = TestFESpace(model, reffe_u)
P = TestFESpace(model, reffe_p)
J = TestFESpace(model, reffe_j)

X = MultiFieldFESpace([U,U])

Ω = Triangulation(model)
dΩ = Measure(Ω, 2)

Γ = Boundary(model)
dΓ = Measure(Γ, 2)

a1(u,v) = ∫(∇(u)⊙∇(v))dΩ
A1 = assemble_matrix(a1,U,U);
@benchmark assemble_matrix!($a1,$A1,$U,$U)

a2(u,v) = ∫(∇(u)⊙∇(v))dΩ + ∫(u⋅v)dΓ
A2 = assemble_matrix(a2,U,U);
@benchmark assemble_matrix!($a2,$A2,$U,$U)

a3(u,v) = ∫(divergence(u)⊙divergence(v))dΩ
A3 = assemble_matrix(a3,U,U);
@benchmark assemble_matrix!($a3,$A3,$U,$U)

a4(u,v) = ∫(divergence(u)⊙divergence(v))dΩ
A4 = assemble_matrix(a4,J,J);
@benchmark assemble_matrix!($a4,$A4,$J,$J)

a5((u1,u2),(v1,v2)) = ∫(∇(u1)⊙∇(v1))dΩ + ∫(∇(u2)⊙∇(v2))dΩ
A5 = assemble_matrix(a5,X,X);
@benchmark assemble_matrix!($a5,$A5,$X,$X)

ProfileView.@profview @benchmark assemble_matrix!($a1,$A1,$U,$U)
ProfileView.@profview @benchmark assemble_matrix!($a2,$A2,$U,$U)
ProfileView.@profview @benchmark assemble_matrix!($a3,$A3,$U,$U)
ProfileView.@profview @benchmark assemble_matrix!($a4,$A4,$J,$J)
