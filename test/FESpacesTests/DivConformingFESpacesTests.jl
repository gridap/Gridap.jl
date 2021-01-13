module DivConformingFESpacesTests

using Test
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Io

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 1

u(x) = x

reffe = ReferenceFE(raviart_thomas,order)

V = TestFESpace(model,reffe,dirichlet_tags = [1,6])
test_single_field_fe_space(V)

U = TrialFESpace(V,u)

uh = interpolate(u,U)

e = u - uh

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

#using Gridap.Visualization
#
#writevtk(Ω,"trian",nsubcells=10,cellfields=["uh"=>uh])

order = 1

reffe = ReferenceFE(TET,raviart_thomas,order)

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = simplexify(CartesianDiscreteModel(domain,partition))

labels = get_face_labeling(model)
dir_tags = Array{Integer}(undef,0)

V = FESpace(model,reffe,conformity=DivConformity())

v(x) = VectorValue(-0.5*x[1]+1.0,-0.5*x[2],-0.5*x[3])
vh = interpolate(v,V)

e = v - vh

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

#using Gridap.Visualization
#
#writevtk(trian,"test",order=3,cellfields=["vh"=>vh])

end # module
