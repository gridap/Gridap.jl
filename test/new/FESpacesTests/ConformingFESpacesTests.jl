module ConformingFESpacesTests

using Test
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry

include("../../../src/new/FESpaces/FESpaces.jl")

using .FESpaces

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

face_labeing = get_face_labeling(model)
dirichlet_tags = ["tag_1","tag_6"]

cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
  reffes, grid_topology, face_labeing, dirichlet_tags)

r = [
  [-1,1,4,5,14,15,16,17,35],[1,2,5,6,18,19,17,20,36],[2,3,6,7,21,22,20,23,37],
  [4,5,8,9,15,24,25,26,38],[5,6,9,10,19,27,26,28,39],[6,7,10,11,22,29,28,30,40],
  [8,9,12,-2,24,-4,31,32,41],[9,10,-2,-3,27,-5,32,33,42],[10,11,-3,13,29,-6,33,34,43]]
test_array(cell_dofs,r)
@test nfree == 43
@test ndiri == 6
@test dirichlet_dof_tag == [1, 2, 2, 2, 2, 2]
@test dirichlet_cells == [1, 7, 8, 9]

end  # module
