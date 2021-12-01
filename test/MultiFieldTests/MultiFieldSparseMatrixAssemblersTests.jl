module MultiFieldSparseMatrixAssemblersTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
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

dy = get_fe_basis(Y)
dx = get_trial_fe_basis(X)
dv, dq = dy
du, dp = dx

grad_dv = ∇(dv)
grad_dp = ∇(dp)

cellmat = integrate(dv*du,quad)
cellvec = integrate(dv*2,quad)
cellmatvec = pair_arrays(cellmat,cellvec)
rows = get_cell_dof_ids(Y,trian)
cols = get_cell_dof_ids(X,trian)
cellmat_c = attach_constraints_cols(X,cellmat,trian)
cellmat_rc = attach_constraints_rows(Y,cellmat_c,trian)
cellvec_r = attach_constraints_rows(Y,cellvec,trian)
cellmatvec_c = attach_constraints_cols(X,cellmatvec,trian)
cellmatvec_rc = attach_constraints_rows(Y,cellmatvec_c,trian)

#cellmat_Γ = integrate(  jump(dv)*dp.⁺ + mean(dq)*jump(dp), quad_Γ)
cellmat_Γ = integrate(  jump(dv)*mean(du) + jump(∇(dq))⋅jump(∇(dp)), quad_Γ)

cellvec_Γ = integrate(  jump(dv) + mean(dq),quad_Γ)
cellmatvec_Γ = pair_arrays(cellmat_Γ,cellvec_Γ)
rows_Γ = get_cell_dof_ids(Y,trian_Γ)
cols_Γ = get_cell_dof_ids(X,trian_Γ)
cellmat_Γ_c = attach_constraints_cols(X,cellmat_Γ,trian_Γ)
cellmat_Γ_rc = attach_constraints_rows(Y,cellmat_Γ_c,trian_Γ)
cellvec_Γ_r = attach_constraints_rows(Y,cellvec_Γ,trian_Γ)
cellmatvec_Γ_c = attach_constraints_cols(X,cellmatvec_Γ,trian_Γ)
cellmatvec_Γ_rc = attach_constraints_rows(Y,cellmatvec_Γ_c,trian_Γ)

assem = SparseMatrixAssembler(SparseMatrixCSR{0,Float64,Int},X,Y)

matvecdata = ([cellmatvec,cellmatvec_Γ],[rows,rows_Γ],[cols,cols_Γ])
matdata = ([cellmat,cellmat_Γ],[rows,rows_Γ],[cols,cols_Γ])
vecdata = ([cellvec,cellvec_Γ],[rows,rows_Γ])
data = (matvecdata,matdata,vecdata)

#@show num_free_dofs(X)
#@show num_free_dofs(Y)
#display(cellmat_Γ[end][1,1])
#display(rows_Γ[end][1])

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
