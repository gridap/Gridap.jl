module RefinementRulesTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.Fields


polys = [QUAD,TRI]

for poly in polys
  rr_generic = RefinementRule(poly,2)
  Adaptivity.test_refinement_rule(rr_generic)

  rr_white = Adaptivity.WhiteRefinementRule(poly)
  Adaptivity.test_refinement_rule(rr_white)

  rr_red = Adaptivity.RedRefinementRule(poly)
  Adaptivity.test_refinement_rule(rr_red)

  n_edges = num_faces(poly,1)
  for e in 1:n_edges
    rr_green = Adaptivity.GreenRefinementRule(poly,e)
    Adaptivity.test_refinement_rule(rr_green)
  end
end

end