module FESpacesWithLinearConstraintsTests

using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Test

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])

V = FESpace(
  model=model,valuetype=Float64,reffe=:Lagrangian,order=1,
  conformity=:H1,dirichlet_tags="dirichlet")

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

@test Vc.n_fdofs == 6
@test Vc.n_fmdofs == 4

fmdof_to_val = collect(Float64,1:num_free_dofs(Vc))
dmdof_to_val = -collect(Float64,1:num_dirichlet_dofs(Vc))
vch = FEFunction(Vc,fmdof_to_val,dmdof_to_val)
r = [[-1.0, -1.5, 1.0, 1.0], [-1.5, -2.0, 1.0, 2.0], [1.0, 1.0, 3.0, 3.5], [1.0, 2.0, 3.5, 4.0]]
@test get_cell_values(vch) ≈ r

u(x) = sin(4*x[1]+0.4)*cos(5*x[2]+0.7)
vch = interpolate(Vc,u)

#trian = Triangulation(model)
#using Gridap.Visualization
#writevtk(trian,"trian",nsubcells=10,cellfields=["vh"=>vh,"vch"=>vch])



end # module
