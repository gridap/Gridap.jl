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

@test get_cell_dof_ids(V) == collect1d(get_cell_node_ids(grid))
@test V.metadata.free_dof_to_node == collect(1:num_nodes(grid))
@test V.metadata.free_dof_to_comp == ones(Int,num_nodes(grid))
@test V.metadata.node_and_comp_to_dof == V.metadata.free_dof_to_node

u(x) = x[1]+x[2]
uh = interpolate(u,V)
uhx = get_free_dof_values(uh)
ux = u.(x)
@test uhx ≈ ux

V = CLagrangianFESpace(VectorValue{2,Float64},grid)
test_single_field_fe_space(V,matvecdata,matdata,vecdata)

@test get_cell_dof_ids(V) == [[1,3,7,9,2,4,8,10], [3,5,9,11,4,6,10,12], [7,9,13,15,8,10,14,16], [9,11,15,17,10,12,16,18]]
@test V.metadata.free_dof_to_node == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9]
@test V.metadata.free_dof_to_comp == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
@test V.metadata.node_and_comp_to_dof == VectorValue{2,Int}[(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16),(17,18)]

u(x) = VectorValue(x[1]+x[2],x[2])
uh = interpolate(u,V)
uhx = get_free_dof_values(uh)
ux = u.(x)
@test uhx ≈ collect1d(reinterpret(ux))

model = CartesianDiscreteModel(domain,partition)
grid = get_grid(model)
labels = get_face_labeling(model)
tags = [1,2,4,5,8]
node_to_tag = get_face_tag_index(labels,tags,0)
masks1 = Bool[true,true,false,true,true]
V = CLagrangianFESpace(Float64,grid,Vector{Float64},node_to_tag,masks1)
@test get_cell_dof_ids(V) == [[-1, -2, 1, 2], [-2, -3, 2, -4], [1, 2, 3, 4], [2, -4, 4, 5]]
@test V.metadata.free_dof_to_node == [4, 5, 7, 8, 9]
@test V.metadata.free_dof_to_comp == [1, 1, 1, 1, 1]
@test V.metadata.dirichlet_dof_to_node == [1,2,3,6]
@test V.metadata.dirichlet_dof_to_comp == [1,1,1,1]
@test V.metadata.node_and_comp_to_dof == [-1, -2, -3, 1, 2, -4, 3, 4, 5]

masks2 = [(true,true),(true,false),(false,false),(false,true),(true,true)]
V = CLagrangianFESpace(VectorValue{2,Float64},grid,Vector{Float64},node_to_tag,masks2)
@test get_cell_dof_ids(V) == [[-1, 1, 3, 5, -2, -3, 4, 6], [1, -4, 5, -5, -3, 2, 6, -6], [3, 5, 7, 9, 4, 6, 8, 10], [5, -5, 9, 11, 6, -6, 10, 12]]
@test V.metadata.free_dof_to_node == [2, 3, 4, 4, 5, 5, 7, 7, 8, 8, 9, 9]
@test V.metadata.free_dof_to_comp == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
@test V.metadata.dirichlet_dof_to_node == [1, 1, 2, 3, 6, 6]
@test V.metadata.dirichlet_dof_to_comp == [1, 2, 2, 1, 1, 2]
@test V.metadata.node_and_comp_to_dof == VectorValue{2, Int32}[(-1, -2), (1, -3), (-4, 2), (3, 4), (5, 6), (-5, -6), (7, 8), (9, 10), (11, 12)]

reffe = ReferenceFE(lagrangian,Float64,2)
V = FESpace(model,reffe)
@test V.metadata === nothing

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(model,reffe,conformity=:L2)
@test V.metadata === nothing

# Check that the factory uses clagrangian when possible

V = FESpace(model,reffe)
@test V.metadata.node_and_comp_to_dof == [1, 2, 3, 4, 5, 6, 7, 8, 9]

V = FESpace(model,reffe,dirichlet_tags=tags)
@test V.metadata.node_and_comp_to_dof == [-1, -2, -3, 1, 2, -4, 3, 4, -5]

V = FESpace(model,reffe,dirichlet_tags=tags,dirichlet_masks=masks1)
@test V.metadata.node_and_comp_to_dof == [-1, -2, -3, 1, 2, -4, 3, 4, 5]

reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},1)
V = FESpace(model,reffe)
@test V.metadata.node_and_comp_to_dof == VectorValue{2, Int32}[(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12), (13, 14), (15, 16), (17, 18)]

V = FESpace(model,reffe,dirichlet_tags=tags)
@test V.metadata.node_and_comp_to_dof  == VectorValue{2, Int32}[(-1, -2), (-3, -4), (-5, -6), (1, 2), (3, 4), (-7, -8), (5, 6), (7, 8), (-9, -10)]

V = FESpace(model,reffe,dirichlet_tags=tags,dirichlet_masks=masks2)
@test V.metadata.node_and_comp_to_dof == VectorValue{2, Int32}[(-1, -2), (1, -3), (-4, 2), (3, 4), (5, 6), (-5, -6), (7, 8), (9, 10), (11, 12)]

# Do not use CLagrangian for models with periodic boundary conditions

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition,isperiodic=(false,true))
reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(model,reffe)
@test V.metadata === nothing

end # module
