module GridWithFEMapsTests

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

T = eltype(get_node_coordinates(model)) #T = VectorValue{1,Float64};
order = 1; pol = Polytope(HEX_AXIS...)

reffe = LagrangianRefFE(T,pol,order)
Vₕ = FESpace(model,reffe;conformity=:H1)

cell_type = Fill(1,num_cells(model))
_reffes = [reffe]
reffes = lazy_map(Reindex(_reffes),cell_type)
grid_map = GridWithFEMap(model,cell_type,reffes)
grid_map = GridWithFEMap(model,reffe)


Vₕ = FESpace(model,reffe;conformity=:H1)

using Gridap.Arrays
reffe = LagrangianRefFE(T,pol,order)
cell_type = Fill(1,num_cells(model))
_reffes = [reffe]
lazy_map(Reindex(_reffes),cell_type)

# GridWithFEMap(model,reffe)

@test get_node_coordinates(grid_map) ≈ get_node_coordinates(model)

Vₕ = grid_map.fe_sp
d(x) = 2*x # Analytical solution (for Dirichlet data)
Uₕ = TrialFESpace(Vₕ,d)
dh = interpolate_everywhere(d,Uₕ)
# Project the grid onto the FESpace (if required)
# Some code here, evaluate the cell_map onto the
# FESpace nodes and assemble (free/fixed) DOFs in
# the global DOF vector(s)

grid = get_grid(model)
@test add_mesh_displacement!(grid_map,dh)  ≈ 3.0*get_node_coordinates(grid)
@test update_coordinates!(grid_map,dh) ≈ 2.0*get_node_coordinates(grid)

end #module
