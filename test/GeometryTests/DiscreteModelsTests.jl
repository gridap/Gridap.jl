module DiscreteModelsTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.TensorValues: det
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock
using Gridap.Io

domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

model_1d = DiscreteModel(Polytope{1},model)
test_discrete_model(model_1d)

model_2d = DiscreteModel(Polytope{2},model)
test_discrete_model(model_2d)

model = DiscreteModelMock()
test_discrete_model(model)

grid = get_grid(model)
topo = get_grid_topology(model)
labeling = get_face_labeling(model)

@test num_dims(model) == 2
@test num_cell_dims(model) == 2
@test num_point_dims(model) == 2
@test num_faces(model) == 9 + 13 + 5
@test num_faces(model,0) == 9
@test num_faces(model,1) == 13
@test num_faces(model,2) == 5
@test num_cells(model) == 5
@test num_vertices(model) == 9
@test num_edges(model) == 13
@test num_facets(model) == 13
@test num_nodes(model) == 9

@test get_cell_node_ids(model) == get_cell_node_ids(grid)
@test get_node_coordinates(model) == get_node_coordinates(grid)
@test get_cell_type(model) == get_cell_type(grid)
@test get_reffes(model) == get_reffes(grid)
@test get_face_nodes(model,0) == get_faces(topo,0,0)
@test get_face_nodes(model,1) == get_faces(topo,1,0)
@test get_face_nodes(model,2) == get_faces(topo,2,0)
@test get_face_nodes(model) == get_face_vertices(topo)
@test get_face_own_nodes(model,0) == get_faces(topo,0,0)
@test get_face_own_nodes(model,1) == empty_table(num_faces(model,1))
@test get_face_own_nodes(model,2) == empty_table(num_faces(model,2))
r = vcat(get_face_own_nodes(model,0),get_face_own_nodes(model,1),get_face_own_nodes(model,2))
@test get_face_own_nodes(model) == r
@test get_vertex_node(model) == collect(1:num_vertices(model))
@test get_node_face_owner(model) == collect(1:num_nodes(model))
@test get_reffaces(ReferenceFE{0},model) == [VERTEX1]
@test get_reffaces(ReferenceFE{1},model) == [SEG2]
@test get_reffaces(ReferenceFE{2},model) == [QUAD4, TRI3]
@test get_reffaces(model) == [VERTEX1, SEG2, QUAD4, TRI3]
@test get_face_type(model,0) == get_face_type(topo,0)
@test get_face_type(model,1) == get_face_type(topo,1)
@test get_face_type(model,2) == get_face_type(topo,2)
@test get_face_type(model) == get_face_type(topo)
@test get_reffaces_offsets(model) == [0,1,2]

grid0 = Grid(ReferenceFE{0},model)
grid1 = Grid(ReferenceFE{1},model)
grid2 = Grid(ReferenceFE{2},model)
test_grid(grid0)
test_grid(grid1)
test_grid(grid2)
@test num_dims(grid0) == 0
@test num_dims(grid1) == 1
@test num_dims(grid2) == 2

model = DiscreteModelMock()
dict = to_dict(model)
model2 = from_dict(DiscreteModel,dict)
test_discrete_model(model2)

domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

tmodel = simplexify(model)
test_discrete_model(tmodel)

ptmodel = simplexify(model,positive=true)
test_discrete_model(ptmodel)

Jtk = lazy_map(∇,get_cell_map(ptmodel))
X0k = lazy_map(first,get_cell_coordinates(ptmodel))
Jt0k = lazy_map(evaluate,Jtk,X0k)
@test all(>(0),lazy_map(det,Jt0k))

model2 = DiscreteModel(grid,topo,labeling)
test_discrete_model(model2)

reffes = ReferenceFE(model,lagrangian,Float64,1)
@test isa(reffes,AbstractVector{<:ReferenceFE})

d = mktempdir()

filename = joinpath(d,"model.json")
to_json_file(model2,filename)
model3 = DiscreteModelFromFile(filename)
test_discrete_model(model3)

rm(d,recursive=true)

end # module
