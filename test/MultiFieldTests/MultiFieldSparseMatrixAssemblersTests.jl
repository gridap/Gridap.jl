module MultiFieldSparseMatrixAssemblersTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using SparseArrays
using SparseMatricesCSR
using Gridap.CellData
using Gridap.MultiField
using Gridap.ReferenceFEs
using Gridap.TensorValues

using Test

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

degree = order
trian = get_triangulation(model)
quad = CellQuadrature(trian,degree)

trian_Γ = SkeletonTriangulation(model)
quad_Γ = CellQuadrature(trian_Γ,degree)

V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:H1)
Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,order-1),conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

free_values = rand(num_free_dofs(X))
xh = FEFunction(X,free_values)
uh, ph = xh

dy = get_cell_shapefuns(Y)
dx = get_cell_shapefuns_trial(X)
dv, dq = dy
du, dp = dx

cellmat = integrate(dv*du,quad)
cellvec = integrate(dv*2,quad)
cellids = get_cell_to_bgcell(trian)
cellmatvec = pair_arrays(cellmat,cellvec)
rows = get_cell_dof_ids(V,cellids)
cols = get_cell_dof_ids(U,cellids)
cellmat_c = attach_constraints_cols(U,cellmat,cellids)
cellmat_rc = attach_constraints_rows(V,cellmat_c,cellids)
cellvec_r = attach_constraints_rows(V,cellvec,cellids)
cellmatvec_c = attach_constraints_cols(U,cellmatvec,cellids)
cellmatvec_rc = attach_constraints_rows(V,cellmatvec_c,cellids)

#cellmat_Γ = integrate(  jump(dv)*dp.⁺ + mean(dq)*jump(dp), quad_Γ)
cellmat_Γ = integrate(  jump(dv)*mean(du) + jump(∇(dq))⋅jump(∇(dp)), quad_Γ)

cellvec_Γ = integrate(  jump(dv) + mean(dq),quad_Γ)
cellmatvec_Γ = pair_arrays(cellmat_Γ,cellvec_Γ)
cellids_Γ = get_cell_to_bgcell(trian_Γ)
rows_Γ = get_cell_dof_ids(V,cellids_Γ)
cols_Γ = get_cell_dof_ids(U,cellids_Γ)
cellmat_Γ_c = attach_constraints_cols(U,cellmat_Γ,cellids_Γ)
cellmat_Γ_rc = attach_constraints_rows(V,cellmat_Γ_c,cellids_Γ)
cellvec_Γ_r = attach_constraints_rows(V,cellvec_Γ,cellids_Γ)
cellmatvec_Γ_c = attach_constraints_cols(U,cellmatvec_Γ,cellids_Γ)
cellmatvec_Γ_rc = attach_constraints_rows(V,cellmatvec_Γ_c,cellids_Γ)

assem = SparseMatrixAssembler(SparseMatrixCSR{0,Float64,Int},X,Y)

matvecdata = ([cellmatvec,cellmatvec_Γ],[rows,rows_Γ],[cols,cols_Γ])
matdata = ([cellmat,cellmat_Γ],[rows,rows_Γ],[cols,cols_Γ])
vecdata = ([cellvec,cellvec_Γ],[rows,rows_Γ])
data = (matvecdata,matdata,vecdata)

test_assembler(assem,matdata,vecdata,data)

A = assemble_matrix(assem,matdata)

A = allocate_matrix(assem,matdata)

assem = SparseMatrixAssembler(X,Y)
test_assembler(assem,matdata,vecdata,data)

struct AssemblyStrategyMock <: AssemblyStrategy end
FESpaces.row_map(a::AssemblyStrategyMock,row) = row
FESpaces.col_map(a::AssemblyStrategyMock,col) = col
FESpaces.row_mask(a::AssemblyStrategyMock,row) = true
FESpaces.col_mask(a::AssemblyStrategyMock,col) = true

assem = SparseMatrixAssembler(
  SparseMatrixCSC{Float64,Int},
  Vector{Float64},
  X,
  Y,
  AssemblyStrategyMock())
test_assembler(assem,matdata,vecdata,data)

end # module
