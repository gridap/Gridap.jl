module CurlConformingFESpacesTests

using Test
using LinearAlgebra
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields
using Gridap.ReferenceFEs

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2

u(x) = VectorValue(x[1]*x[1],x[1]*x[1]*x[1],0.0)
# u(x) = x

reffe = ReferenceFE(:Nedelec,order)

@test_broken begin

V = TestFESpace(model,reffe,dirichlet_tags = [21,22])
  # dirichlet_tags = "boundary")
test_single_field_fe_space(V)

U = TrialFESpace(V,u)

uh = interpolate(u,U)

e = u - uh

Ω = Triangulation(model)
dΩ = LebesgueMeasure(Ω,order)

el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
el2 < 1.0e-10
end

# using Gridap.Visualization

# writevtk(trian,"trian",nsubcells=10,cellfields=["uh"=>uh])


end # module
