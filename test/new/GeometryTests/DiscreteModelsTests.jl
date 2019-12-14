module DiscreteModelsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock

model = DiscreteModelMock()
test_discrete_model(model)
@test is_oriented(model) == false

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
@test get_vertex_coordinates(model) == get_node_coordinates(model)
@test get_isboundary_face(model) == Bool[1,1,1,1,0,1,1,1,1,1,0,1,0,1,0,1,0,1,1,0,1,1,1,1,1,1,1]
@test get_isboundary_node(model) == Bool[1,1,1,1,0,1,1,1,1]
@test get_dimranges(model) == [1:9, 10:22, 23:27]
@test get_offsets(model) == [0, 9, 22]
@test get_offset(model,0) == 0
@test get_offset(model,1) == 9
@test get_offset(model,2) == 22
@test get_facedims(model) == Int8[0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2]
@test get_cell_faces(model) == [
  [1,2,4,5,10,11,12,13,23], [2,3,5,14,13,15,24], [3,6,5,16,15,17,25],
  [4,5,7,8,11,18,19,20,26], [5,6,8,9,17,21,20,22,27]]

grid = UnstructuredGrid(ReferenceFE{0},model)
@test num_cells(grid) == num_vertices(model)

grid = ConformingTriangulation(ReferenceFE{1},model)
@test num_cells(grid) == num_edges(model)

grid = ConformingTriangulation(ReferenceFE{2},model)
@test num_cells(grid) == num_cells(model)

# TODO for the moment only implemented for the aligned case
#order = 2
#reffes = [ LagrangianRefFE(Float64,get_polytope(reffe),order) for reffe in get_reffes(model)]
#model2 = replace_reffes(model,reffes)
#test_discrete_model(model2)

end # module
