module FESpacesWithLinearConstraintsTests

using Gridap.Algebra
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.FESpaces
using Test
using LinearAlgebra

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
trian = Triangulation(model)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[6,7,8])

V = FESpace(
  model=model,valuetype=Float64,reffe=:Lagrangian,order=1,
  conformity=:H1,dirichlet_tags="dirichlet")
test_single_field_fe_space(V)

fdof_to_val = collect(Float64,1:num_free_dofs(V))
ddof_to_val = -collect(Float64,1:num_dirichlet_dofs(V))
vh = FEFunction(V,fdof_to_val,ddof_to_val)

sDOF_to_dof = [1,5,-2]
sDOF_to_dofs = Table([[-1,4],[4,6],[-1,-3]])
sDOF_to_coeffs = Table([[0.5,0.5],[0.5,0.5],[0.5,0.5]])

Vc = FESpaceWithLinearConstraints(
  sDOF_to_dof,
  sDOF_to_dofs,
  sDOF_to_coeffs,
  V)

test_single_field_fe_space(Vc)
@test has_constraints(Vc)

@test Vc.n_fdofs == 6
@test Vc.n_fmdofs == 4

fmdof_to_val = collect(Float64,1:num_free_dofs(Vc))
dmdof_to_val = -collect(Float64,1:num_dirichlet_dofs(Vc))
vch = FEFunction(Vc,fmdof_to_val,dmdof_to_val)
r = [[-1.0, -1.5, 1.0, 1.0], [-1.5, -2.0, 1.0, 2.0], [1.0, 1.0, 3.0, 3.5], [1.0, 2.0, 3.5, 4.0]]
@test get_cell_values(vch) ≈ r

v(x) = sin(4*x[1]+0.4)*cos(5*x[2]+0.7)
vch = interpolate(v,Vc)

#using Gridap.Visualization
#writevtk(trian,"trian",nsubcells=10,cellfields=["vh"=>vh,"vch"=>vch])

u(x) = x[1] + 2*x[2]
f(x) = -Δ(u)(x)
Uc = TrialFESpace(Vc,u)
@test has_constraints(Uc)
uch = interpolate(u,Uc)

btrian = BoundaryTriangulation(model,"neumann")

quad = CellQuadrature(trian,2)
bquad = CellQuadrature(btrian,2)

bn = get_normal_vector(btrian)

a(u,v) = ∇(v)⋅∇(u)
b1(v) = v*f
b2(v) = v*(bn⋅∇(u))
t1 = AffineFETerm(a,b1,trian,quad)
t2 = FESource(b2,btrian,bquad)
op = AffineFEOperator(Uc,Vc,t1,t2)
uch = solve(op)

#using Gridap.Visualization
#writevtk(trian,"trian",nsubcells=10,cellfields=["uch"=>uch])

e = u - uch

l2(u) = u*u
h1(u) = u*u + a(u,u)

e_l2 = sqrt(sum(integrate(l2(e),trian,quad)))
e_h1 = sqrt(sum(integrate(h1(e),trian,quad)))

tol = 1.e-9
@test e_l2 < tol
@test e_h1 < tol

end # module
