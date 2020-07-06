module CDLagrangianFESpacesTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
trian = get_triangulation(model)

# order = 1
orders = (2,1)
order = 2*max(orders...)

T = VectorValue{2,Float64}
reffe = CDLagrangianRefFE(T,QUAD,orders,(DISC,CONT))
V = FESpace(model=model,reffe=reffe)#),dirichlet_tags = [1,6])
test_single_field_fe_space(V)

u(x) = x

U = TrialFESpace(V,u)
uh = interpolate(U,u)

e = u - uh

trian = Triangulation(model)
quad = CellQuadrature(trian,order)

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

reffe = CDLagrangianRefFE(Float64,QUAD,(2,2),(CONT,DISC))

quad9 = LagrangianRefFE(T,QUAD,2)
V = FESpace(model=model,reffe=quad9,conformity=CDConformity((CONT,DISC)))

U = TrialFESpace(V,u)
uh = interpolate(U,u)
e = u - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

cdquad9 = CDLagrangianRefFE(T,QUAD,2,(CONT,DISC))
V = FESpace(model=model,reffe=cdquad9)

U = TrialFESpace(V,u)
uh = interpolate(U,u)
e = u - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

cdquad9 = CDLagrangianRefFE(T,QUAD,2,(CONT,DISC))
V = FESpace(model=model,reffe=cdquad9,conformity=CDConformity((CONT,DISC)))

U = TrialFESpace(V,u)
uh = interpolate(U,u)
e = u - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

cdquad9 = CDLagrangianRefFE(T,QUAD,2,(CONT,DISC))
V = FESpace(model=model,reffe=cdquad9,conformity=CDConformity((DISC,CONT)))

U = TrialFESpace(V,u)
uh = interpolate(U,u)
e = u - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

cdquad9 = LagrangianRefFE(T,QUAD,2,(CONT,DISC))
V = FESpace(model=model,reffe=cdquad9,conformity=CDConformity((DISC,CONT)))

U = TrialFESpace(V,u)
uh = interpolate(U,u)
e = u - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

cdquad9 = LagrangianRefFE(T,QUAD,(2,1))
V = FESpace(model=model,reffe=cdquad9,conformity=CDConformity((DISC,CONT)))

U = TrialFESpace(V,u)
uh = interpolate(U,u)
e = u - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

cdquad9 = LagrangianRefFE(T,QUAD,(2,0))
V = FESpace(model=model,reffe=cdquad9,conformity=CDConformity((CONT,DISC)))

u(x) = VectorValue(x[1],0.0)
U = TrialFESpace(V,u)
uh = interpolate(U,u)
e = u - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
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
