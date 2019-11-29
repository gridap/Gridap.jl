module DiscreteModelsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: ConformingTrianMock
using Gridap.Geometry: DiscreteModelMock

import Gridap.ReferenceFEs
import Gridap.ReferenceFEs: num_faces
import Gridap.ReferenceFEs: num_nodes
import Gridap.ReferenceFEs: get_faces
import Gridap.ReferenceFEs: get_face_nodes

import Gridap.Geometry: get_vertex_node
import Gridap.Geometry: get_node_face_owner
import Gridap.Geometry: get_face_nodes
import Gridap.Geometry: get_cell_nodes
import Gridap.Geometry: get_isboundary_face
import Gridap.Geometry: get_face_reffe_type
import Gridap.Geometry: get_face_polytope_type
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_polytopes
import Gridap.Geometry: get_node_coordinates

using Gridap.Helpers
using Gridap.Arrays
using Gridap.ReferenceFEs: _find_unique_with_indices
include("../../../src/new/Geometry/GridOperations.jl")

include("../../../src/new/Geometry/UnstructuredDiscreteModels.jl")

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

end # module
