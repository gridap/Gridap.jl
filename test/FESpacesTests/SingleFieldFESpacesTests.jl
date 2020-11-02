module SingleFieldFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
reffes = ReferenceFE(model,basis=:Lagrangian,order=order,valuetype=Float64)
V0 = FESpace(model,reffes,dirichlet_tags=["tag_24","tag_25"])

matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(V0,matvecdata,matdata,vecdata)

f(x) = sin(4*pi*(x[1]-x[2]^2))+1

fh = interpolate_everywhere(f, V0)

fh = interpolate_dirichlet(f, V0)

dirichlet_values = compute_dirichlet_values_for_tags(V0,[1,2])

r=[2,2,1,1,2,2,1,1,2,2,2,2,1,2,1,1,1,1,2,2,2,2,1,2,1,1,1,1,2,2,1,1,2,2,1,2,1,1,2,2,1,2,1,1,2,2,1,2,1,1]

@test dirichlet_values == r

dirichlet_values = compute_dirichlet_values_for_tags(V0,2)
@test dirichlet_values == fill(2,num_dirichlet_dofs(V0))

free_values = zero_free_values(V0)
uh = FEFunction(V0,free_values,dirichlet_values)


#using Gridap.Visualization
#
#trian = get_triangulation(model)
#
#writevtk(model,"model")
#
#writevtk(trian,"trian",cellfields=["fh" => fh, "uh" => uh])

#
#
#strian = SkeletonTriangulation(model)
#
#fh_Γ = restrict(get_array(fh),strian)
#
#writevtk(strian,"strian",nsubcells=10,cellfields=["jump_fh" => jump(fh_Γ)])

end # module
