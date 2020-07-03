module FESpaceFactoriesTests

using Test

using Gridap.Geometry
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.FESpaces

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
order = 3

T = VectorValue{2,Float64}
reffe = QDiscRefFE(T,QUAD,2)
V = FESpace(model=model,reffe=reffe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model=model,reffe=QUAD4,conformity=GradConformity())
v1 = FEFunction(V,rand(num_free_dofs(V)))
@test isa(V,UnconstrainedFESpace)

V = FESpace(model=model,reffe=QUAD4,conformity=L2Conformity())
v2 = FEFunction(V,rand(num_free_dofs(V)))
@test isa(V,UnconstrainedFESpace)

V = FESpace(model=model,reffe=QUAD4,conformity=CDConformity((CONT,DISC)))
v3 = FEFunction(V,rand(num_free_dofs(V)))
@test isa(V,UnconstrainedFESpace)

using Gridap.Visualization
writevtk(Triangulation(model),"results",cellfields=["v1"=>v1,"v2"=>v2,"v3"=>v3])

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:Lagrangian,
 order=order,
 conformity=:L2)

@test isa(V,UnconstrainedFESpace)

V = FESpace(
 triangulation=get_triangulation(model),
 valuetype=Float64,
 reffe=:Lagrangian,
 order=order,
 conformity=:L2)

@test isa(V,UnconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:SLagrangian,
 order=order,
 conformity=:H1)

@test isa(V,UnconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:QLagrangian,
 conformity=:H1,
 order=order,
 dirichlet_tags="boundary")

@test isa(V,UnconstrainedFESpace)

D = num_point_dims(model)

V = FESpace(
 model=model,
 valuetype=VectorValue{D,Float64},
 reffe=:Lagrangian,
 order=order,
 conformity=:H1,
 dirichlet_tags="boundary",
 dirichlet_masks=(true,false))

@test isa(V,UnconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=VectorValue{D,Float64},
 reffe=:SLagrangian,
 order=order,
 conformity=:H1,
 dirichlet_tags=[1,2],
 dirichlet_masks=[(true,false),(false,true)])

@test isa(V,UnconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:PLagrangian,
 order=order,
 conformity=:L2)

#using Gridap.Visualization
#uh = FEFunction(V,rand(num_free_dofs(V)))
#writevtk(get_triangulation(model),"trian",nsubcells=20,cellfields=["uh"=>uh])

@test isa(V,UnconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:PLagrangian,
 order=order,
 conformity=:L2,
 constraint=:zeromean)

@test isa(V,ZeroMeanFESpace)

trian = reindex(Triangulation(model),1:4)
quad = CellQuadrature(trian,order)

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:PLagrangian,
 order=order,
 conformity=:L2,
 constraint=:zeromean,
 zeromean_trian = trian,
 zeromean_quad = quad)

@test abs(V.vol - 4/9) < 1.0e-9

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:PLagrangian,
 order=order,
 conformity=:L2,
 constraint=:zeromean,
 zeromean_trian = trian)

@test abs(V.vol - 4/9) < 1.0e-9

end # module
