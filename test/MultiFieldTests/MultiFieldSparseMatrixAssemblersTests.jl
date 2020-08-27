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

V = TestFESpace(model=model,order=order,reffe=:Lagrangian,conformity=:H1,valuetype=Float64)
Q = TestFESpace(model=model,order=order-1,reffe=:Lagrangian,conformity=:L2,valuetype=Float64)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

free_values = rand(num_free_dofs(X))
xh = FEFunction(X,free_values)
uh, ph = xh

dy = get_cell_basis(Y)
dx = get_cell_basis(X)
dv, dq = dy
du, dp = dx

cellmat = integrate(dv*du,trian,quad)
cellvec = integrate(dv*2,trian,quad)
cellids = get_cell_id(trian)
cellmatvec = pair_arrays(cellmat,cellvec)

dv_Γ, dq_Γ = restrict(get_cell_basis(Y), trian_Γ)
du_Γ, dp_Γ = restrict(get_cell_basis(X), trian_Γ)
cellmat_Γ = integrate(  jump(dv_Γ)*dp_Γ.⁺ + mean(dq_Γ)*jump(dp_Γ), trian_Γ,quad_Γ)
cellvec_Γ = integrate(  jump(dv_Γ) + mean(dq_Γ), trian_Γ,quad_Γ)
cellmatvec_Γ = pair_arrays(cellmat_Γ,cellvec_Γ)
cellids_Γ = get_cell_id(trian_Γ)

assem = SparseMatrixAssembler(SparseMatrixCSR{0},X,Y)

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
