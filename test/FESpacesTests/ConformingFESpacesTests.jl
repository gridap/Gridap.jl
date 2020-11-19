module ConformingFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields

# testing compute_conforming_cell_dofs

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]
cell_reffe = expand_cell_data(reffes,get_cell_type(grid_topology))

face_labeling = get_face_labeling(model)
dirichlet_tags = ["tag_1","tag_6"]

trian = Triangulation(model)
cell_map = get_cell_map(trian)
cell_fe = CellFE(cell_map,cell_reffe)

cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
  cell_fe,CellConformity(cell_fe),grid_topology, face_labeling, dirichlet_tags)

r = [
  [-1,1,4,5,14,15,16,17,35],[1,2,5,6,18,19,17,20,36],[2,3,6,7,21,22,20,23,37],
  [4,5,8,9,15,24,25,26,38],[5,6,9,10,19,27,26,28,39],[6,7,10,11,22,29,28,30,40],
  [8,9,12,-2,24,-4,31,32,41],[9,10,-2,-3,27,-5,32,33,42],[10,11,-3,13,29,-6,33,34,43]]
test_array(cell_dofs,r)
@test nfree == 43
@test ndiri == 6
@test dirichlet_dof_tag == [1, 2, 2, 2, 2, 2]
@test dirichlet_cells == [1, 7, 8, 9]

order = 1
reffes = [LagrangianRefFE(VectorValue{2,Float64},p,order) for p in polytopes]
cell_reffe = expand_cell_data(reffes,get_cell_type(grid_topology))

dirichlet_components = [(true,true), (false,true)]

cell_fe = CellFE(cell_map,cell_reffe)

cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
  cell_fe,CellConformity(cell_fe),grid_topology, face_labeling, dirichlet_tags, dirichlet_components)

r = [
  [-1,1,7,9,-2,2,8,10],[1,3,9,11,2,4,10,12],[3,5,11,13,4,6,12,14],
  [7,9,15,17,8,10,16,18],[9,11,17,19,10,12,18,20],[11,13,19,21,12,14,20,22],
  [15,17,23,25,16,18,24,-3],[17,19,25,26,18,20,-3,-4],[19,21,26,27,20,22,-4,28]]

test_array(cell_dofs,r)
@test nfree==28
@test ndiri==4
@test dirichlet_dof_tag == [1, 1, 2, 2,]
@test dirichlet_cells == [1, 7, 8, 9]

order = 3
reffes = [LagrangianRefFE(VectorValue{2,Float64},p,order) for p in polytopes]
cell_reffe = expand_cell_data(reffes,get_cell_type(grid_topology))

dirichlet_components = [(true,true), (false,true)]
cell_fe = CellFE(cell_map,cell_reffe)

cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
  cell_fe,CellConformity(cell_fe), grid_topology, face_labeling, dirichlet_tags, dirichlet_components)

reffe = ReferenceFE(:Lagrangian,VectorValue{2,Float64},3)

V = FESpace(model,reffe,dirichlet_tags=dirichlet_tags)

test_single_field_fe_space(V)

matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(V,matvecdata,matdata,vecdata)

V = FESpace(model,reffe,dirichlet_tags=dirichlet_tags,dirichlet_masks=dirichlet_components)
test_single_field_fe_space(V)

end  # module
