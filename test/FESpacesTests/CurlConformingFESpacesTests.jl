module CurlConformingFESpacesTests

using Test
using LinearAlgebra
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields
using Gridap.ReferenceFEs

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
order = 1
u((x,y)) = 2*VectorValue(2,3)
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10
#using Gridap.Visualization
#writevtk(Ω,"nedel",nsubcells=10,cellfields=["err"=>e,"u"=>u,"uh"=>uh])

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition) |> simplexify
order = 0
u((x,y)) = 2*VectorValue(-y,x)
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

order = 1
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

order = 2
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition) |> simplexify
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
u((x,y,z)) = 2*VectorValue(-y,x,0.) - VectorValue(0.,-z,y)
order = 0
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10
order = 1
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10
order = 2
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2

u(x) = VectorValue(x[1]*x[1],x[1]*x[1]*x[1],0.0)
#u(x) = VectorValue(2,3,5)
# u(x) = x

reffe = ReferenceFE(nedelec,order)


V = TestFESpace(model,reffe,dirichlet_tags = [21,22])
# dirichlet_tags = "boundary")
test_single_field_fe_space(V)

U = TrialFESpace(V,u)

uh = interpolate(u,U)

e = u - uh

Ω = Triangulation(model)
dΩ = Measure(Ω,order)

el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10


# using Gridap.Visualization

# writevtk(trian,"trian",nsubcells=10,cellfields=["uh"=>uh])


end # module
