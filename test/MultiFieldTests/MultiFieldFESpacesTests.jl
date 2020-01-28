module MultiFieldFESpacesTests

using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using Test

using Gridap.MultiField
using Gridap.MultiField: MultiFieldFESpace
using Gridap.MultiField: ConsequtiveMultiFieldStyle

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

multi_field_style = ConsequtiveMultiFieldStyle()

Y = MultiFieldFESpace([V,Q],multi_field_style)
X = MultiFieldFESpace([U,P],multi_field_style)

@test num_free_dofs(X) == num_free_dofs(U) + num_free_dofs(P)
@test num_free_dofs(X) == num_free_dofs(Y)

free_values = rand(num_free_dofs(X))
xh = FEFunction(X,free_values)
@test is_a_fe_function(xh)
uh, ph = xh
@test is_a_fe_function(uh)
@test is_a_fe_function(ph)

dy = get_cell_basis(Y)
@test is_test(dy)
@test is_a_fe_cell_basis(dy)
dv, dq = dy
@test is_a_fe_cell_basis(dv)
@test is_a_fe_cell_basis(dq)

dx = get_cell_basis(X)
@test is_a_fe_cell_basis(dx)
@test is_trial(dx)
du, dp = dx
@test is_a_fe_cell_basis(du)
@test is_a_fe_cell_basis(dp)

cellmat = integrate(dv*du,trian,quad)

cellvec = integrate(dv*2,trian,quad)

cellids = get_cell_id(trian)

test_fe_space(V,cellmat,cellvec,cellids,cellids)
test_fe_space(U,cellmat,cellvec,cellids,cellids)

#using Gridap.Visualization
#writevtk(trian,"trian";nsubcells=30,cellfields=["uh" => uh, "ph"=> ph])


end # module
