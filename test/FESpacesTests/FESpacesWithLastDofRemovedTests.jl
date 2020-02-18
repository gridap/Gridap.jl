module FESpacesWithLastDofRemovedTests

using Gridap.Geometry
using Gridap.FESpaces

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)

order = 2

V = TestFESpace(
  model=model,
  reffe=:Lagrangian,
  valuetype=Float64,
  order=order,
  conformity=:L2)

V0 = FESpaceWithLastDofRemoved(V)
test_single_field_fe_space(V0)

fun(x) = sin(4*pi*(x[1]+x[2]^2)) + 3
uh0 = interpolate(V0,fun)

#using Gridap.Visualization
#
#writevtk(trian,"trian",nsubcells=20,cellfields=["uh0"=>uh0])

end # module
