module FESpaceFactoriesTests

using Test

using Gridap.Geometry
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.FESpaces

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
order = 3

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:Lagrangian,
 order=order,
 conformity=:L2)

@test isa(V,UnsconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:SLagrangian,
 order=order,
 conformity=:H1)

@test isa(V,UnsconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=Float64,
 reffe=:QLagrangian,
 conformity=:H1,
 order=order,
 dirichlet_tags="boundary")

@test isa(V,UnsconstrainedFESpace)

D = num_point_dims(model)

V = FESpace(
 model=model,
 valuetype=VectorValue{D,Float64},
 reffe=:Lagrangian,
 order=order,
 conformity=:H1,
 dirichlet_tags="boundary",
 dirichlet_masks=(true,false))

@test isa(V,UnsconstrainedFESpace)

V = FESpace(
 model=model,
 valuetype=VectorValue{D,Float64},
 reffe=:SLagrangian,
 order=order,
 conformity=:H1,
 dirichlet_tags=[1,2],
 dirichlet_masks=[(true,false),(false,true)])

@test isa(V,UnsconstrainedFESpace)

end # module
