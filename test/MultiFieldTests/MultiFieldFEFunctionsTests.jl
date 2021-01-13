module MultiFieldFEFunctionsTests

using BlockArrays
using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using Gridap.CellData
using Gridap.MultiField
using Gridap.ReferenceFEs
using Test

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:H1)
Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,order-1),conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

free_values = rand(num_free_dofs(X))
xh = FEFunction(X,free_values)
test_fe_function(xh)
uh, ph = xh

cell_values = get_cell_dof_values(xh)
@test isa(cell_values[1],BlockArrayCoo)

end # module
