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
reffe = ReferenceFE(:Lagrangian,order=order,valuetype=Float64)
V = FESpace(model,reffe,dirichlet_tags=["tag_1","tag_6"])
test_fe_space(V)

vh = FEFunction(V,rand(num_free_dofs(V)))
@test isa(vh,FEFunction)
test_fe_function(vh)

cellids = [1,3,5,2]
cell_vals = get_cell_dof_values(vh,cellids)
@test cell_vals == get_cell_dof_values(vh)[cellids]

cellidsL = cellids
cellidsR = [2,4,3,1]
cellidsS = SkeletonPair(cellidsL,cellidsR)
@test_broken begin
cell_vals = get_cell_dof_values(vh,cellidsS)
isa(cell_vals[1],BlockArrayCoo)
end

zh = zero(V)
@test isa(zh,FEFunction)
@test has_constraints(V) == false

cellids = [1,3,5,2]
@test get_cell_dof_ids(V,cellids) == [[-1, 1, 4, 5], [2, 3, 6, 7], [5, 6, 9, 10], [1, 2, 5, 6]]

cellidsL = cellids
cellidsR = [2,4,3,1]
cellidsS = SkeletonPair(cellidsL,cellidsR)
@test_broken isa(get_cell_dof_ids(V,cellidsS)[1],BlockArrayCoo)
@test_broken get_cell_dof_ids(V,cellidsS)[1] == [-1, 1, 4, 5, 1, 2, 5, 6]

cell_constr = get_cell_constraints(V)
@test cell_constr == [Matrix(I,4,4) for cell in 1:num_cells(model)]

cell_constr = get_cell_constraints(V,cellids)
@test cell_constr == [Matrix(I,4,4) for cell in 1:length(cellids)]

@test get_cell_isconstrained(V,cellids) == Fill(false,length(cellids))

@test_broken begin
cell_constr = get_cell_constraints(V,cellidsS)
isa(cell_constr[1],BlockArrayCoo)
end

du = get_cell_shapefuns_trial(V)
du_data = get_cell_data(du)
@test size(du_data[1]) == (1,4)

@test get_cell_isconstrained(V,cellidsS) == Fill(false,length(cellids))

cellmat = [rand(4,4) for cell in 1:num_cells(model)]
cellvec = [rand(4) for cell in 1:num_cells(model)]
cellmatvec = pair_arrays(cellmat,cellvec)
cellids = IdentityVector(num_cells(model))
matdata = (cellmat,cellids,cellids)
vecdata = (cellvec,cellids)
matvecdata = (cellmatvec,cellids,cellids)

@test cellmat === attach_constraints_rows(V,cellmat,cellids)
@test cellmat === attach_constraints_cols(V,cellmat,cellids)
test_fe_space(V,matvecdata,matdata,vecdata)

end # module
