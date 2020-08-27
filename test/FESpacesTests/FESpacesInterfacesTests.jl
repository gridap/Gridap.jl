module FESpacesInterfacesTests

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
test_fe_space(V)

zh = zero(V)
@test is_a_fe_function(zh)
@test has_constraints(V) == false

cellids = [1,3,5,2]
@test get_cell_dofs(V,cellids) == [[-1, 1, 4, 5], [2, 3, 6, 7], [5, 6, 9, 10], [1, 2, 5, 6]]

cellidsL = cellids
cellidsR = [2,4,3,1]
cellidsS = SkeletonPair(cellidsL,cellidsR)
@test isa(get_cell_dofs(V,cellidsS)[1],BlockArrayCoo)
@test get_cell_dofs(V,cellidsS)[1] == [-1, 1, 4, 5, 1, 2, 5, 6]

@test get_cell_axes(V) == get_cell_axes_with_constraints(V)

cell_constr = get_cell_constraints(V)
@test cell_constr == [Matrix(I,4,4) for cell in 1:num_cells(model)]

cell_constr = get_cell_constraints(V,cellids)
@test cell_constr == [Matrix(I,4,4) for cell in 1:length(cellids)]

cell_constr = get_cell_constraints(V,cellidsS)
@test isa(cell_constr[1],BlockArrayCoo)

@test get_cell_isconstrained(V,cellids) == Fill(false,length(cellids))

@test get_cell_isconstrained(V,cellidsS) == Fill(false,length(cellids))

cellmat = [rand(4,4) for cell in 1:num_cells(model)]
cellvec = [rand(4) for cell in 1:num_cells(model)]
cellmatvec = pair_arrays(cellmat,cellvec)
cellids = identity_vector(num_cells(model))
matdata = (cellmat,cellids,cellids)
vecdata = (cellvec,cellids)
matvecdata = (cellmatvec,cellids,cellids)

@test cellmat === attach_constraints_rows(V,cellmat,cellids)
@test cellmat === attach_constraints_cols(V,cellmat,cellids)
test_fe_space(V,matvecdata,matdata,vecdata)

end # module
