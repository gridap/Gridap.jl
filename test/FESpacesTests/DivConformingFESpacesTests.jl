module DivConformingFESpacesTests

using Test
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Fields

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
trian = get_triangulation(model)

order = 1

u(x) = x

V = TestFESpace(
  reffe = :RaviartThomas,
  conformity = :Hdiv,
  order = order,
  model = model,
  dirichlet_tags = [1,6])
test_single_field_fe_space(V)

U = TrialFESpace(V,u)

uh = interpolate(u,U)

e = u - uh

trian = Triangulation(model)
quad = CellQuadrature(trian,order)

el2 = sqrt(sum(integrate(inner(e,e),quad)))
@test el2 < 1.0e-10

T = Float64
reffe = RaviartThomasRefFE(T,QUAD,order)
V = FESpace(model=model,reffe=reffe,dirichlet_tags = [1,6])
test_single_field_fe_space(V)

U = TrialFESpace(V,u)

uh = interpolate(u,U)

el2 = sqrt(sum(integrate(inner(e,e),quad)))
@test el2 < 1.0e-10

#using Gridap.Visualization
#
#writevtk(trian,"trian",nsubcells=10,cellfields=["uh"=>uh])

end # module
