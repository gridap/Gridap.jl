module MultiFieldSparseMatrixAssemblersTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using SparseArrays
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

#cellmat_Γ = integrate(  jump(dv)*dp.⁺ + mean(dq)*jump(dp), quad_Γ)
cellmat_Γ = integrate(  jump(dv)*mean(du) + jump(∇(dq))⋅jump(∇(dp)), quad_Γ)

cellvec_Γ = integrate(  jump(dv) + mean(dq),quad_Γ)
cellmatvec_Γ = pair_arrays(cellmat_Γ,cellvec_Γ)
cellids_Γ = get_cell_to_bgcell(trian_Γ)

assem = SparseMatrixAssembler(SparseMatrixCSR{0,Float64,Int},X,Y)

matvecdata = ([cellmatvec,cellmatvec_Γ],[cellids,cellids_Γ],[cellids,cellids_Γ])
matdata = ([cellmat,cellmat_Γ],[cellids,cellids_Γ],[cellids,cellids_Γ])
vecdata = ([cellvec,cellvec_Γ],[cellids,cellids_Γ])
data = (matvecdata,matdata,vecdata)

test_assembler(assem,matdata,vecdata,data)

A = assemble_matrix(assem,matdata)

A = allocate_matrix(assem,matdata)

assem = SparseMatrixAssembler(X,Y)
test_assembler(assem,matdata,vecdata,data)

end # module
