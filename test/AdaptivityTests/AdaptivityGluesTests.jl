module AdaptivityGluesTests

using Test
using Gridap
using Gridap.Io
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

rr     = RefinementRule(QUAD, (2,2))
rrules = Fill(rr, 4)
glue   = Adaptivity.blocked_refinement_glue(rrules)

# in-memory dict round-trip
d     = to_dict(glue)
glue2 = from_dict(AdaptivityGlue, d)
@test glue.n2o_faces_map[end]   == glue2.n2o_faces_map[end]
@test glue.n2o_cell_to_child_id == glue2.n2o_cell_to_child_id
@test collect(glue.o2n_faces_map) == collect(glue2.o2n_faces_map)

# JSON round-trip
glue3 = from_json(AdaptivityGlue, to_json(glue))
@test glue.n2o_faces_map[end]   == glue3.n2o_faces_map[end]
@test glue.n2o_cell_to_child_id == glue3.n2o_cell_to_child_id

end # module
