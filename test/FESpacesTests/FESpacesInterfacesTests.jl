module FESpacesInterfacesTests

using FillArrays
using Test
using Gridap.Arrays
using Gridap.Fields
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
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model,reffe,dirichlet_tags=["tag_1","tag_6"])
test_fe_space(V)


vh = FEFunction(V,rand(num_free_dofs(V)))
@test isa(vh,FEFunction)
test_fe_function(vh)


dv = get_fe_basis(V)
du = get_trial_fe_basis(V)

trian_Γ = SkeletonTriangulation(model)
x_Γ = get_cell_points(trian_Γ)

@test isa(dv.minus(x_Γ)[1],ArrayBlock)
@test isa(du.plus(x_Γ)[1],ArrayBlock)
@test isa(∇(dv).plus(x_Γ)[1],ArrayBlock)
@test isa(∇(du).minus(x_Γ)[1],ArrayBlock)

cellids = [1,3,5,2]
cell_vals = get_cell_dof_values(vh,cellids)
@test cell_vals == get_cell_dof_values(vh)[cellids]

cellidsL = cellids
cellidsR = [2,4,3,1]
cellidsS = SkeletonPair(cellidsL,cellidsR)
cell_vals = get_cell_dof_values(vh,cellidsS)
@test isa(cell_vals[1],ArrayBlock)

zh = zero(V)
@test isa(zh,FEFunction)
@test has_constraints(V) == false

cellids = [1,3,5,2]
@test get_cell_dof_ids(V,cellids) == [[-1, 1, 4, 5], [2, 3, 6, 7], [5, 6, 9, 10], [1, 2, 5, 6]]

cellidsL = cellids
cellidsR = [2,4,3,1]
cellidsS = SkeletonPair(cellidsL,cellidsR)
@test isa(get_cell_dof_ids(V,cellidsS)[1],ArrayBlock)
@test get_cell_dof_ids(V,cellidsS)[1][1] == [-1, 1, 4, 5]
@test get_cell_dof_ids(V,cellidsS)[1][2] == [1, 2, 5, 6]

cell_constr = get_cell_constraints(V)
@test cell_constr == [Matrix(I,4,4) for cell in 1:num_cells(model)]

cell_constr = get_cell_constraints(V,cellids)
@test cell_constr == [Matrix(I,4,4) for cell in 1:length(cellids)]

@test get_cell_isconstrained(V,cellids) == Fill(false,length(cellids))

cell_constr = get_cell_constraints(V,cellidsS)
@test isa(cell_constr[1],ArrayBlock)

du = get_trial_fe_basis(V)
du_data = get_data(du)
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
