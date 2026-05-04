module HighOrderGridsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields

# Helpers

function _num_ho_nodes(partition, order, ::Val{D}) where D
  # For an n-cube grid with uniform order, the number of high-order nodes
  # equals the number of nodes of a refined grid with (order*partition) cells
  # at order 1, i.e. prod(order .* partition .+ 1).
  prod(order .* partition .+ 1)
end

# ---- 2D quad, order 2 ----

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain, partition)
grid  = UnstructuredGrid(get_grid(model))
topo  = get_grid_topology(model)

ho = high_order_grid(grid, topo, 2)
@test num_cells(ho) == prod(partition)
@test num_nodes(ho) == _num_ho_nodes(partition, 2, Val(2))

# Node coordinates must lie inside [0,1]^2
coords = get_node_coordinates(ho)
@test all(p -> all(xi -> 0 <= xi <= 1, Tuple(p)), coords)

# ---- 2D simplex, order 3 ----

model_s = simplexify(CartesianDiscreteModel(domain, partition))
grid_s  = UnstructuredGrid(get_grid(model_s))
topo_s  = get_grid_topology(model_s)

ho_s = high_order_grid(grid_s, topo_s, 3)
@test num_cells(ho_s) == num_cells(grid_s)
coords_s = get_node_coordinates(ho_s)
@test all(p -> all(xi -> -1e-14 <= xi <= 1+1e-14, Tuple(p)), coords_s)

# ---- 3D hex, order 2 ----

domain3 = (0,1,0,1,0,1)
partition3 = (2,2,2)
model3 = CartesianDiscreteModel(domain3, partition3)
grid3  = UnstructuredGrid(get_grid(model3))
topo3  = get_grid_topology(model3)

ho3 = high_order_grid(grid3, topo3, 2)
@test num_cells(ho3) == prod(partition3)
@test num_nodes(ho3) == _num_ho_nodes(partition3, 2, Val(3))
coords3 = get_node_coordinates(ho3)
@test all(p -> all(xi -> 0 <= xi <= 1, Tuple(p)), coords3)

# ---- Order-1 roundtrip ----
# A first-order high_order_grid should recover the same number of nodes
# and (up to reordering) the same coordinates.

ho1 = high_order_grid(grid, topo, 1)
@test num_nodes(ho1) == num_nodes(grid)
@test num_cells(ho1) == num_cells(grid)

orig = sort(get_node_coordinates(grid), by=Tuple)
ho1c = sort(get_node_coordinates(ho1),  by=Tuple)
@test all(map((a,b) -> isapprox(a, b; atol=1e-14), orig, ho1c))

# ---- Explicit cell_maps keyword ----
# Passing cell_maps explicitly should give the same result as the default.

cell_maps = get_cell_map(grid)
ho_explicit = high_order_grid(grid, topo, 2; cell_maps)
@test num_nodes(ho_explicit) == num_nodes(ho)
coords_expl = get_node_coordinates(ho_explicit)
coords_ho   = get_node_coordinates(ho)
@test all(map((a,b) -> isapprox(a, b; atol=1e-14), coords_expl, coords_ho))

end # module
