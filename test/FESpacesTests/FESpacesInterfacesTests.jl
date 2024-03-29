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

Ω = Triangulation(model)
∂Ω = BoundaryTriangulation(model)
trian_Γ = SkeletonTriangulation(model)
x_Γ = get_cell_points(trian_Γ)

@test isa(dv.minus(x_Γ)[1],ArrayBlock)
@test isa(du.plus(x_Γ)[1],ArrayBlock)
@test isa(∇(dv).plus(x_Γ)[1],ArrayBlock)
@test isa(∇(du).minus(x_Γ)[1],ArrayBlock)

glue = get_glue(∂Ω,Val(num_cell_dims(model)))
cellids = glue.tface_to_mface
cell_vals = get_cell_dof_values(vh,∂Ω)
@test cell_vals == get_cell_dof_values(vh)[cellids]

cell_vals = get_cell_dof_values(vh,trian_Γ)
@test isa(cell_vals[1],ArrayBlock)

zh = zero(V)
@test isa(zh,FEFunction)
@test has_constraints(V) == false

cellids = [1,3,5,2]
Ω1 = Triangulation(model,cellids)
@test get_cell_dof_ids(V,Ω1) == [[-1, 1, 4, 5], [2, 3, 6, 7], [5, 6, 9, 10], [1, 2, 5, 6]]

@test isa(get_cell_dof_ids(V,trian_Γ)[1],ArrayBlock)
@test get_cell_dof_ids(V,trian_Γ)[1][1] == [-1, 1, 4, 5]
@test get_cell_dof_ids(V,trian_Γ)[1][2] == [4, 5, 8, 9]

cell_constr = get_cell_constraints(V)
@test cell_constr == [Matrix(I,4,4) for cell in 1:num_cells(model)]

cell_constr = get_cell_constraints(V,Ω1)
@test cell_constr == [Matrix(I,4,4) for cell in 1:length(cellids)]

@test get_cell_isconstrained(V,Ω1) == Fill(false,length(cellids))

cell_constr = get_cell_constraints(V,trian_Γ)
@test isa(cell_constr[1],ArrayBlock)

du = get_trial_fe_basis(V)
du_data = get_data(du)
@test size(du_data[1]) == (1,4)

@test get_cell_isconstrained(V,trian_Γ) == Fill(false,num_cells(trian_Γ))

cellmat = [rand(4,4) for cell in 1:num_cells(model)]
cellvec = [rand(4) for cell in 1:num_cells(model)]
cellmatvec = pair_arrays(cellmat,cellvec)

@test cellmat === attach_constraints_rows(V,cellmat,Ω)
@test cellmat === attach_constraints_cols(V,cellmat,Ω)
test_fe_space(V,cellmatvec,cellmat,cellvec,Ω)

end # module
