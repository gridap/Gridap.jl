module RefinementRuleBoundaryTests

using Test
using Gridap
using Gridap.Helpers
using Gridap.Adaptivity
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays

# Test the final function in 2D
rr = Gridap.Adaptivity.RedRefinementRule(QUAD)
res = Adaptivity.get_face_subface_ldof_to_cell_ldof(rr,(2,2),1)

# Test the final function in 3D
rr = Gridap.Adaptivity.RedRefinementRule(HEX)
res = Adaptivity.get_face_subface_ldof_to_cell_ldof(rr,(2,2,2),1)
res = Adaptivity.get_face_subface_ldof_to_cell_ldof(rr,(2,2,2),2)

end