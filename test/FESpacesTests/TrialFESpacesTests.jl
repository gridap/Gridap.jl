module TrialFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.FESpaces

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

dirichlet_tags = ["tag_01","tag_10"]
V = GradConformingFESpace(reffes,model,dirichlet_tags)

U = TrialFESpace(V,[4,3])
test_single_field_fe_space(U,[],[],[],[])

uh = interpolate(U,0)

@test get_dirichlet_values(U) == [4.0, 3.0, 3.0, 3.0, 3.0, 3.0]

cell_basis = get_cell_basis(U)
@test is_trial(cell_basis)

#trian = get_triangulation(model)
#
#using Gridap.Visualization
#
#writevtk(trian,"trian",cellfields=["uh"=>uh])

end # module
