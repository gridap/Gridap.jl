module DiscreteModelWithFEMapsTests

using Test
using Gridap
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.CellData
using FillArrays

# Create a model

domain = (0,1)
partition = (5,)
model = CartesianDiscreteModel(domain,partition)

order = 1;
grid_map = GridWithFEMap(model,order)

orders = Fill(1,num_cells(model))
grid_map = GridWithFEMap(model,orders)

@test get_node_coordinates(grid_map) ≈ get_node_coordinates(model)

T = eltype(get_node_coordinates(model))
pol = Polytope(HEX_AXIS...)
reffe = LagrangianRefFE(T,pol,order)
Vₕ = FESpace(model,reffe;conformity=:H1)

Vₕ = grid_map.fe_sp
d(x) = 2*x # Analytical solution (for Dirichlet data)
Uₕ = TrialFESpace(Vₕ,d)
dh = interpolate_everywhere(d,Uₕ)

grid = get_grid(model)

add_mesh_displacement!(grid_map,dh)
@test get_node_coordinates(grid_map)  ≈ 3.0*get_node_coordinates(grid)

update_coordinates!(grid_map,dh)
@test get_node_coordinates(grid_map) ≈ 2.0*get_node_coordinates(grid)

end #module
