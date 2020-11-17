module CDLagrangianFESpacesTests

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

# order = 1
orders = (2,1)
order = 2*max(orders...)

Ω = Triangulation(model)
dΩ = LebesgueMeasure(Ω,order)

T = VectorValue{2,Float64}

V = FESpace(model,ReferenceFE(:Lagrangian,T,orders),conformity=CDConformity((DISC,CONT)))
test_single_field_fe_space(V)
u(x) = x
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum(∫(e⋅e)*dΩ))
@test el2 < 1.0e-10

reffe = LagrangianRefFE(T,QUAD,2)
V = FESpace(model,ReferenceFE(:Lagrangian,T,2),conformity=CDConformity((CONT,DISC)))
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum(∫(e⋅e)*dΩ))
@test el2 < 1.0e-10

V = FESpace(model,ReferenceFE(:Lagrangian,T,(2,0)),conformity=CDConformity((CONT,DISC)))
u(x) = VectorValue(x[1],0.0)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum(∫(e⋅e)*dΩ))
@test el2 < 1.0e-10

V = FESpace(model,ReferenceFE(:Lagrangian,T,(2,0)),conformity=CDConformity((DISC,DISC)))
u(x) = VectorValue(x[1],0.0)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum(∫(e⋅e)*dΩ))
@test el2 < 1.0e-10


# # Skeleton triangulation
# face_own_dofs = get_face_own_dofs(reffe)
# strian = SkeletonTriangulation(model,reffe,face_own_dofs)
# ns = get_normal_vector(strian)
# writevtk(strian,"strian",cellfields=["normal"=>ns])
#
# # Random function for visualization purposes
# model = CartesianDiscreteModel((0,1,0,1),(10,5))
# V = FESpace(model=model,reffe=reffe,conformity=true)
# trian = Triangulation(model)
# vh = FEFunction(V,rand(num_free_dofs(V)))
# writevtk(trian,"trian",nsubcells=20,cellfields=["vh"=>vh])

#using Gridap.Visualization
#
#writevtk(trian,"trian",nsubcells=10,cellfields=["uh"=>uh])

end # module
