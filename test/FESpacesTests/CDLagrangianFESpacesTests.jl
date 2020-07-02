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

#using Gridap.Visualization
#
#writevtk(trian,"trian",nsubcells=10,cellfields=["uh"=>uh])

end # module
