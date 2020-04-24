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

fdof_to_dofs =   Table([[-1,4],[2],[3],[4],[4,6],[6]])
fdof_to_coeffs = Table([[0.5,0.5],[1.],[1.],[1.],[0.5,0.5],[1.]])
ddof_to_dofs = Table([[-1],[-1,-3],[-3]])
ddof_to_coeffs = Table([[1.],[0.5,0.5],[1.]])

Vc = FESpaceWithLinearConstraints(
  fdof_to_dofs,
  fdof_to_coeffs,
  ddof_to_dofs,
  ddof_to_coeffs,
  V)

@test Vc.n_fdofs == 6
@test Vc.n_fmdofs == 4

fmdof_to_val = collect(Float64,1:num_free_dofs(Vc))
dmdof_to_val = -collect(Float64,1:num_dirichlet_dofs(Vc))
vch = FEFunction(Vc,fmdof_to_val,dmdof_to_val)
u(x) = sin(4*x[1]+0.4)*cos(5*x[2]+0.7)
vch = interpolate(Vc,u)

display(get_cell_values(vch))

trian = Triangulation(model)
using Gridap.Visualization
writevtk(trian,"trian",nsubcells=10,cellfields=["vh"=>vh,"vch"=>vch])

#display(Vc.mDOF_to_DOF)
#display(Vc.DOF_to_mDOFs)
#display(Vc.DOF_to_coeffs)
#display(Vc.cell_to_lmdof_to_mdof)


end # module
