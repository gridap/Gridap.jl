module MultiFieldSparseMatrixAssemblersTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using SparseArrays

using Test

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

V = TestFESpace(model=model,order=order,reffe=:Lagrangian,conformity=:H1,valuetype=Float64)
Q = TestFESpace(model=model,order=order-1,reffe=:Lagrangian,conformity=:L2,valuetype=Float64)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = [V,Q]
X = [U,P]

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

assem = SparseMatrixAssembler(SparseMatrixCSR{0},Y,X)
test_assembler(assem,[cellmat],[cellvec],[cellids],[cellids])

A = assemble_matrix(assem,[cellmat],[cellids],[cellids])

A = allocate_matrix(assem,[cellids],[cellids])

assem = SparseMatrixAssembler(Y,X)
test_assembler(assem,[cellmat],[cellvec],[cellids],[cellids])

end # module
