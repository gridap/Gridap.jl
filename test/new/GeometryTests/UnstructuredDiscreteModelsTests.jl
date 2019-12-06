module UnstructuredDiscreteModelsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: ConformingTrianMock
using Gridap.Geometry: DiscreteModelMock

grid = ConformingTrianMock()

model = UnstructuredDiscreteModel(grid)
test_discrete_model(model)

m = DiscreteModelMock()

@test num_faces(model,0) == num_faces(m,0)
@test num_faces(model,1) == num_faces(m,1)
@test num_faces(model,2) == num_faces(m,2)
@test num_nodes(model) == num_nodes(m)
@test get_reffes(ReferenceFE{1},model) == get_reffes(ReferenceFE{1},m)
@test get_reffes(ReferenceFE{0},model) == get_reffes(ReferenceFE{0},m)
@test get_reffes(ReferenceFE{2},model) == get_reffes(ReferenceFE{2},m)
@test get_polytopes(Polytope{1},model) == get_polytopes(Polytope{1},m)
@test get_polytopes(Polytope{0},model) == get_polytopes(Polytope{0},m)
@test get_polytopes(Polytope{2},model) == get_polytopes(Polytope{2},m)
@test get_face_reffe_type(model,1) == get_face_reffe_type(m,1)
@test get_face_polytope_type(model,1) == get_face_polytope_type(m,1)
@test get_isboundary_face(model,0) == get_isboundary_face(m,0)
@test get_isboundary_face(model,1) == get_isboundary_face(m,1)
@test get_isboundary_face(model,2) == get_isboundary_face(m,2)
@test get_face_nodes(model,1) == get_face_nodes(m,1)
@test get_face_nodes(model,0) == get_face_nodes(m,0)
@test get_face_nodes(model,2) == get_face_nodes(m,2)
@test get_isboundary_node(model) == get_isboundary_node(m)
@test get_faces(model,2,0) == get_faces(m,2,0)
@test get_faces(model,2,1) == get_faces(m,2,1)
@test get_faces(model,1,2) == get_faces(m,1,2)
@test get_faces(model,2,1) == get_faces(m,2,1)
@test get_faces(model,0,1) == get_faces(m,0,1)
@test get_faces(model,1,0) == get_faces(m,1,0)
@test get_faces(model,1,1) == get_faces(m,1,1)
@test get_node_coordinates(model) == get_node_coordinates(m)
@test get_vertex_coordinates(model) == get_vertex_coordinates(m)           

labels = get_face_labeling(model)
@test num_entities(labels) == 0
@test num_tags(labels) == 0
@test get_face_labeling(model) === labels
@test get_face_entity(get_face_labeling(model),0) === get_face_entity(labels,0)

domain = (0,1,0,1,0,1)
partition = (2,2,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)

edge_to_isboundary = get_isboundary_face(model,1)
@test length(edge_to_isboundary) == num_edges(model)

edge_to_faces = get_faces(model,1,2)
@test length(edge_to_faces) == num_edges(model)
@test maximum(edge_to_faces.data) == num_facets(model)
face_to_edges = get_faces(model,2,1)
@test length(face_to_edges) == num_facets(model)
@test maximum(face_to_edges.data) == num_edges(model)

test_discrete_model(model)

domain = (0,1,0,1,0,1)
partition = (2,3,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)
test_discrete_model(model)

end # module
