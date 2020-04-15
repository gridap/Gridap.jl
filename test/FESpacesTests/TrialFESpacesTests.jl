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
@test U.dirichlet_values == compute_dirichlet_values_for_tags(V,[4,3])
v = copy(U.dirichlet_values)
v .= 0
isa(V,SingleFieldFESpace)
ud = compute_dirichlet_values_for_tags!(v,V,[4,3])
@test all(ud .== v)
test_single_field_fe_space(U)
U = TrialFESpace!(v,V,[4,3])


matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(U,matvecdata,matdata,vecdata)

uh = interpolate(U,0)

@test get_dirichlet_values(U) == [4.0, 3.0, 3.0, 3.0, 3.0, 3.0]
TrialFESpace!(U,[1,2])
@test get_dirichlet_values(U) == [1.0, 2.0, 2.0, 2.0, 2.0, 2.0]

cell_basis = get_cell_basis(U)
@test is_trial(cell_basis)

U0 = HomogeneousTrialFESpace(V)
@test get_dirichlet_values(U0) == zeros(6)

U0 = HomogeneousTrialFESpace!(v,V)
@test v === get_dirichlet_values(U0)
@test v == zeros(6)
@test get_dirichlet_values(U0) == zeros(6)

#trian = get_triangulation(model)
#
#using Gridap.Visualization
#
#writevtk(trian,"trian",cellfields=["uh"=>uh])

end # module
