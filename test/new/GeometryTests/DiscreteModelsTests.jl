module DiscreteModelsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry

using Gridap.Helpers

import Gridap.Geometry: num_cell_dims
import Gridap.Geometry: num_point_dims

import Gridap.ReferenceFEs: get_node_coordinates
import Gridap.ReferenceFEs: num_nodes
import Gridap.ReferenceFEs: is_affine
import Gridap.ReferenceFEs: has_straight_faces
import Gridap.ReferenceFEs: get_reffes
import Gridap.ReferenceFEs: get_faces
import Gridap.ReferenceFEs: num_dims
import Gridap.ReferenceFEs: get_dimranges
import Gridap.ReferenceFEs: num_faces
import Gridap.ReferenceFEs: num_facets
import Gridap.ReferenceFEs: num_edges
import Gridap.ReferenceFEs: num_vertices
import Gridap.ReferenceFEs: get_facedims
import Gridap.ReferenceFEs: get_offsets
import Gridap.ReferenceFEs: get_offset
import Gridap.ReferenceFEs: get_vertex_coordinates

include("../../../src/new/Geometry/DiscreteModels.jl")

include("../../../src/new/Geometry/DiscreteModelMocks.jl")

model = DiscreteModelMock()
test_discrete_model(model)

@test num_dims(model) == 2
@test num_cell_dims(model) == 2
@test num_point_dims(model) == 2

end # module
