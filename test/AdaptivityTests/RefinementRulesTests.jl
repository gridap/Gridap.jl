module RefinementRulesTests

using Test
using Gridap
using Gridap.Adaptivity
using Gridap.ReferenceFEs

polys = [TRI,QUAD]

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

rr_bc2 = Adaptivity.BarycentricRefinementRule(TRI)
Adaptivity.test_refinement_rule(rr_bc2)

rr_bc3 = Adaptivity.BarycentricRefinementRule(TET)
Adaptivity.test_refinement_rule(rr_bc3)

rr_ps2 = Adaptivity.PowellSabinRefinementRule(TRI)
rr_ps3 = Adaptivity.PowellSabinRefinementRule(TET)

end