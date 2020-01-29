module DiscontinuousFESpacesTests

using Test
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 3
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

trian = get_triangulation(model)

V = DiscontinuousFESpace(reffes,trian)
test_single_field_fe_space(V,[],[],[],[])

U = TrialFESpace(V)

f(x) = sin(pi*x[1])*cos(2*pi*x[2])

fh = interpolate(U,f)

uh = FEFunction(V,rand(num_free_dofs(V)))

#using Gridap.Visualization
#
#writevtk(trian,"trian",nsubcells=40,cellfields=["fh"=>fh, "uh"=>uh])

reffes = [PDiscRefFE(Float64,p,order) for p in polytopes]
V = DiscontinuousFESpace(reffes,trian)
test_single_field_fe_space(V,[],[],[],[])

end # module
