module CurlConformingFESpacesTests
##
using Test
using Gridap
using Gridap.MapApply
using Gridap.CellValuesGallery
using LinearAlgebra

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
trian = get_triangulation(model)

order = 2

u(x) = x

V = TestFESpace(
  reffe = :Nedelec,
  conformity = :Hcurl,
  order = order,
  model = model,
  dirichlet_tags = [1,6])
test_single_field_fe_space(V)

U = TrialFESpace(V,u)

uh = interpolate(U,u)

# #using Gridap.Visualization
# #
# #writevtk(trian,"trian",nsubcells=10,cellfields=["uh"=>uh])
# D = 2
# order = 2
#
# dom = fill(1,D)
# part = Tuple(fill(5,D))
# model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=part)
# # model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=part)
#
# p = Polytope(fill(HEX_AXIS,D)...)
# reffe = NedelecRefFE(p,Float64,order)
#
#
# grid = Grid(model,D)
# trian = Triangulation(grid)
# graph = GridGraph(model)
#
#
#
# labels = FaceLabels(model)
# # tags = [5,6,7,8]
# tags = Int[]
# fun(x) = VectorValue(x[2],x[1],x[2]*x[1])
# # fun(x) = VectorValue(1.0,1.0)
# # tags = Int[]
#
# fesp = ConformingFESpace(reffe,trian,graph,labels,tags)
# dfesp = TrialFESpace(fesp,fun)
#
#
# trian = Triangulation(model)
# quad = CellQuadrature(trian,degree=2)
#
# uh = interpolate(dfesp,fun)
# uh.free_dofs
# uh.diri_dofs
#
# e = fun - uh
# uh.free_dofs
# uh.diri_dofs
#
#
# el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
# @test el2 < 1.0e-10

##

end # module
