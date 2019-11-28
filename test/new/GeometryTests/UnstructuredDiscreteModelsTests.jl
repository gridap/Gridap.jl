module DiscreteModelsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: ConformingTrianMock

using Gridap.Helpers
using Gridap.Arrays

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

using Gridap.ReferenceFEs: _find_unique_with_indices
include("../../../src/new/Geometry/GridOperations.jl")

include("../../../src/new/Geometry/UnstructuredDiscreteModels.jl")

grid = ConformingTrianMock()

model = UnstructuredDiscreteModel(grid)

@show num_faces(model,0)
@show num_faces(model,1)
@show num_faces(model,2)
@show num_nodes(model)
@show get_reffes(ReferenceFE{1},model)
@show get_face_reffe_type(model,1)
@show get_isboundary_face(model,1)

#@show model
#
@show get_faces(model,2,0)
@show get_faces(model,2,1)
@show get_faces(model,1,2)
@show get_faces(model,2,1)
@show get_faces(model,0,1)
@show get_faces(model,1,0)

@show model

end # module
