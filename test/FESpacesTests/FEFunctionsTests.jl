module FEFunctionsTests

using FillArrays
using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using LinearAlgebra

order = 1
domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]
dirichlet_tags = ["tag_1","tag_6"]
V = GradConformingFESpace(reffes,model,dirichlet_tags)

vh = FEFunction(V,rand(num_free_dofs(V)))
@test is_a_fe_function(vh)
test_fe_function(vh)

cellids = [1,3,5,2]
cell_vals = get_cell_values(vh,cellids)
@test cell_vals == reindex(get_cell_values(vh),cellids)

cellidsL = cellids
cellidsR = [2,4,3,1]
cellidsS = SkeletonPair(cellidsL,cellidsR)
cell_vals = get_cell_values(vh,cellidsS)
@test isa(cell_vals[1],BlockArrayCoo)

end # module
