module MultiFieldFESpacesWithLinearConstraintsTests

using Gridap.Algebra
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.TensorValues
using Gridap.MultiField
using Gridap.CellData
using Gridap.ReferenceFEs
using Test

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_1",[1,2,5])
add_tag_from_tags!(labels,"neumann_1",[6,7,8])
add_tag_from_tags!(labels,"dirichlet_2",[1,3,7])
add_tag_from_tags!(labels,"neumann_2",[5,6,8])

Ω = Triangulation(model)
Γ1 = BoundaryTriangulation(model,tags="neumann_1")
Γ2 = BoundaryTriangulation(model,tags="neumann_2")
Λ = SkeletonTriangulation(model)

dΩ = LebesgueMeasure(Ω,2)
dΓ1 = LebesgueMeasure(Γ1,2)
dΓ2 = LebesgueMeasure(Γ2,2)
dΛ = LebesgueMeasure(Λ,2)

n_Γ1 = get_normal_vector(Γ1)
n_Γ2 = get_normal_vector(Γ2)

reffe = ReferenceFE(:Lagrangian,Float64,1)
V1 = FESpace(model,reffe,conformity=:H1,dirichlet_tags="dirichlet_1")
V2 = FESpace(model,reffe,conformity=:H1,dirichlet_tags="dirichlet_2")

sDOF_to_dof = [1,5,-2]
sDOF_to_dofs = Table([[-1,4],[4,6],[-1,-3]])
sDOF_to_coeffs = Table([[0.5,0.5],[0.5,0.5],[0.5,0.5]])
V1c = FESpaceWithLinearConstraints(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V1)

sDOF_to_dof = [4,-2]
sDOF_to_dofs = Table([[2,6],[-1,-3]])
sDOF_to_coeffs = Table([[0.5,0.5],[0.5,0.5]])
V2c = FESpaceWithLinearConstraints(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V2)

u1(x) = x[1] + 2*x[2]
u2(x) = 4*x[1] + x[2]
f1(x) = -Δ(u1)(x) + u2(x)
#f2(x) = -Δ(u2)(x)

U1c = TrialFESpace(V1c,u1)
U2c = TrialFESpace(V2c,u2)

U = MultiFieldFESpace([U1c,U2c])
V = MultiFieldFESpace([V1c,V2c])
test_fe_space(U)
test_fe_space(V)
@test has_constraints(U)
@test has_constraints(V)

a((u1,u2),(v1,v2)) =
  ∫( ∇(v1)⋅∇(u1) + ∇(v2)⋅∇(u2) + v1*u2 )*dΩ +
  ∫( jump(v1)*jump(u2) + jump(v2)*jump(u2) )*dΛ

b((v1,v2)) =
  ∫( v1*f1 )*dΩ +
  ∫( v1*(n_Γ1⋅∇(u1)) )*dΓ1 +
  ∫( v2*(n_Γ2⋅∇(u2)) )*dΓ2 +
  ∫( jump(v2) )*dΛ

op = AffineFEOperator(a,b,U,V)
uh = solve(op)

uh1, uh2 = uh
e1 = u1 - uh1
e2 = u2 - uh2

#using Gridap.Visualization
#v1 = FEFunction(V1,rand(num_free_dofs(V1)))
#v2 = FEFunction(V2,rand(num_free_dofs(V2)))
#v1c = FEFunction(V1c,rand(num_free_dofs(V1c)))
#v2c = FEFunction(V2c,rand(num_free_dofs(V2c)))
#writevtk(Ω,"trian",nsubcells=10,cellfields=["v1"=>v1,"v2"=>v2,"v1c"=>v1c,"v2c"=>v2c])

#using Gridap.Visualization
#writevtk(Ω,"trian",nsubcells=10,cellfields=["uh1"=>uh1,"uh2"=>uh2,"e1"=>e1,"e2"=>e2])

e1_l2 = sqrt(sum(∫(e1*e1)*dΩ))
e2_l2 = sqrt(sum(∫(e2*e2)*dΩ))

tol = 1.0e-9
@test e1_l2 < tol
@test e2_l2 < tol


end # module
