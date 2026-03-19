module RefinementRulesTests

using Test
using Gridap
using Gridap.Io
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

# IO round-trip
rr_io = RefinementRule(QUAD, (2,2))

d   = to_dict(rr_io)
rr2 = from_dict(RefinementRule, d)
@test get_polytope(rr_io) == get_polytope(rr2)
@test Adaptivity.num_subcells(rr_io) == Adaptivity.num_subcells(rr2)
@test Adaptivity.RefinementRuleType(rr_io) == Adaptivity.RefinementRuleType(rr2)

rr3 = from_json(RefinementRule, to_json(rr_io))
@test get_polytope(rr_io) == get_polytope(rr3)
@test Adaptivity.num_subcells(rr_io) == Adaptivity.num_subcells(rr3)

end