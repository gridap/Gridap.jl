module FESpacesWithLinearConstraintsTests

using Test
using LinearAlgebra
using Gridap

using Gridap.Algebra
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.ReferenceFEs

model = CartesianDiscreteModel((0,1,0,1),(2,2))

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[6,7,8])

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model,tags="neumann")
Λ = SkeletonTriangulation(model)
dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
dΛ = Measure(Λ,2)

V = FESpace(
  model,ReferenceFE(lagrangian,Float64,1);
  conformity=:H1, dirichlet_tags="dirichlet"
)
test_single_field_fe_space(V)

fdof_to_val = collect(Float64,1:num_free_dofs(V))
ddof_to_val = -collect(Float64,1:num_dirichlet_dofs(V))
vh = FEFunction(V,fdof_to_val,ddof_to_val)

sDOF_to_dof = [1,5,-2]
sDOF_to_dofs = Table([[-1,4],[4,6],[-1,-3]])
sDOF_to_coeffs = Table([[0.5,0.5],[0.5,0.5],[0.5,0.5]])

Vc = FESpaceWithLinearConstraints(
  sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V
)

test_single_field_fe_space(Vc)
@test has_constraints(Vc)

@test isa(get_cell_constraints(Vc,Λ)[1],ArrayBlock)
@test num_free_dofs(Vc) == 4

fmdof_to_val = collect(Float64,1:num_free_dofs(Vc))
dmdof_to_val = -collect(Float64,1:num_dirichlet_dofs(Vc))
vch = FEFunction(Vc,fmdof_to_val,dmdof_to_val)
r = [[-1.0, -1.5, 1.0, 1.0], [-1.5, -2.0, 1.0, 2.0], [1.0, 1.0, 3.0, 3.5], [1.0, 2.0, 3.5, 4.0]]
@test get_cell_dof_values(vch) ≈ r

v(x) = sin(4*x[1]+0.4)*cos(5*x[2]+0.7)
vch = interpolate(v,Vc)

#using Gridap.Visualization
#writevtk(Ω,"trian",nsubcells=10,cellfields=["vh"=>vh,"vch"=>vch])

u(x) = x[1] + 2*x[2]
f(x) = -Δ(u)(x)
Uc = TrialFESpace(Vc,u)
@test has_constraints(Uc)
uch = interpolate(u,Uc)

n_Γ = get_normal_vector(Γ)
a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ + ∫( jump(u)*jump(v) )*dΛ
l(v) = ∫( v*f )*dΩ + ∫( v*(n_Γ⋅∇(u)) )*dΓ

op = AffineFEOperator(a,l,Uc,Vc)
uch = solve(op)

#using Gridap.Visualization
#writevtk(trian,"trian",nsubcells=10,cellfields=["uch"=>uch])

e = u - uch
e_l2 = sqrt(sum(∫( e*e )*dΩ))
e_h1 = sqrt(sum(∫( e*e + ∇(e)⋅∇(e) )*dΩ))

tol = 1.e-9
@test e_l2 < tol
@test e_h1 < tol

# Test with complex values

V2 = FESpace(
  model,ReferenceFE(lagrangian,Float64,1);
  conformity=:H1, dirichlet_tags="dirichlet", vector_type=ComplexF64
)

Vc2 = FESpaceWithLinearConstraints(
  sDOF_to_dof,
  sDOF_to_dofs,
  sDOF_to_coeffs,
  V2
)

@test get_dof_value_type(Vc2) <: ComplexF64

# Alternative constructor

fdof_to_dofs = Table([[-1,4],[2],[3],[4],[4,6],[6]])
fdof_to_coeffs = Table([[0.5,0.5],[1.0],[1.0],[1.0],[0.5,0.5],[1.0]])
ddof_to_dofs = Table([[-1],[-1,-3],[-3]])
ddof_to_coeffs = Table([[1.0],[0.5,0.5],[1.0]])

Vc3 = FESpaceWithLinearConstraints(
  fdof_to_dofs, fdof_to_coeffs, ddof_to_dofs, ddof_to_coeffs, V
)

test_single_field_fe_space(Vc3)

@test Vc3.mDOF_to_dof == Vc.mDOF_to_dof
@test Vc3.sDOF_to_dof == Vc.sDOF_to_dof
@test Vc3.sDOF_to_mdofs == Vc.sDOF_to_mdofs
@test Vc3.sDOF_to_coeffs == Vc.sDOF_to_coeffs

end # module
