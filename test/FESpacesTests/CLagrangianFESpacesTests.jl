module CLagrangianFESpacesTests

using Test
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.CellData

domain = (0,1,0,1)
partition = (2,2)
grid = CartesianGrid(domain,partition)
x = collect1d(get_node_coordinates(grid))

V = CLagrangianFESpace(Float64,grid)
matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(V,matvecdata,matdata,vecdata)

@test get_cell_dof_ids(V) === get_cell_node_ids(grid)
@test V.dof_to_node == collect(1:num_nodes(grid))
@test V.dof_to_comp == ones(Int,num_nodes(grid))
@test V.node_and_comp_to_dof == V.dof_to_node

u(x) = x[1]+x[2]
uh = interpolate(u,V)
uhx = get_free_values(uh)
ux = u.(x)
@test uhx ≈ ux

V = CLagrangianFESpace(VectorValue{2,Float64},grid)
test_single_field_fe_space(V,matvecdata,matdata,vecdata)

@test get_cell_dof_ids(V) == [[1,3,7,9,2,4,8,10], [3,5,9,11,4,6,10,12], [7,9,13,15,8,10,14,16], [9,11,15,17,10,12,16,18]]
@test V.dof_to_node == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9]
@test V.dof_to_comp == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
@test V.node_and_comp_to_dof == VectorValue{2,Int}[(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16),(17,18)]

u(x) = VectorValue(x[1]+x[2],x[2])
uh = interpolate(u,V)
uhx = get_free_values(uh)
ux = u.(x)
@test uhx ≈ collect1d(reinterpret(ux))

end # module
