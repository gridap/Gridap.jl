module ZeroMeanFESpacesTests

using Test
using Gridap.Geometry
using Gridap.Fields
using Gridap.Integration
using Gridap.FESpaces
using Gridap.CellData

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 2

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)


_V = TestFESpace(
  model=model,
  reffe=:Lagrangian,
  valuetype=Float64,
  order=order,
  conformity=:L2)

V = ZeroMeanFESpace(_V,trian,quad)
matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(V,matvecdata,matdata,vecdata)

U = TrialFESpace(V)
test_single_field_fe_space(U,matvecdata,matdata,vecdata)
@test isa(U,ZeroMeanFESpace)
@test is_trial(get_cell_basis(U))

fun(x) = sin(4*pi*(x[1]+x[2]^2)) + 3
uh = interpolate(fun, U)

mean1 = sum(integrate(uh,quad))

tol = 1.0e-10
@test abs(mean1) < tol

end # module
