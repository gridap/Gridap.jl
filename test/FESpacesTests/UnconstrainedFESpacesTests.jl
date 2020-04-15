module UnconstrainedFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces

domain =(0,1,0,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

dirichlet_tags = ["tag_1","tag_6"]
V0 = GradConformingFESpace(reffes,model,dirichlet_tags)

f(x) = sin(4*pi*(x[1]-x[2]^2))+1

fh = interpolate(V0,f)
@test get_dirichlet_values(fh) == get_dirichlet_values(V0)

fh = interpolate_everywhere(V0,f)

fh = interpolate_dirichlet(V0,f)
@test get_free_values(fh) == zero_free_values(V0)

#trian = get_triangulation(model)
#
#using Gridap.Visualization
#
#writevtk(trian,"trian",nsubcells=10,cellfields=["fh" => fh])

end # module
