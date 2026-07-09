module AdaptivityCompositionTests

using Test
using Gridap
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays

rr_white = Adaptivity.WhiteRefinementRule(TRI)
rr_red   = Adaptivity.RedRefinementRule(TRI)
rr_bary  = Adaptivity.BarycentricRefinementRule(TRI)

############################################################################################
# compose_refinement_rules: White identities
############################################################################################

@test Adaptivity.compose_refinement_rules(rr_white, rr_red)  === rr_red
@test Adaptivity.compose_refinement_rules(rr_red,   rr_white) === rr_red
@test Adaptivity.compose_refinement_rules(rr_white, rr_bary) === rr_bary
@test Adaptivity.compose_refinement_rules(rr_bary,  rr_white) === rr_bary
@test Adaptivity.RefinementRuleType(Adaptivity.compose_refinement_rules(rr_white, rr_white)) isa Adaptivity.WithoutRefinement

############################################################################################
# compose_refinement_rules: Red ∘ Red
############################################################################################

rr_red2 = Adaptivity.compose_refinement_rules(rr_red, rr_red)
@test Adaptivity.RefinementRuleType(rr_red2) isa Adaptivity.GenericRefinement
@test Adaptivity.get_polytope(rr_red2) == TRI
@test Adaptivity.num_subcells(rr_red2) == Adaptivity.num_subcells(rr_red)^2

############################################################################################
# compose_refinement_rules: generic fallback (Bary ∘ Red and Red ∘ Bary)
############################################################################################

# Bary ∘ Red: Red splits into 4, Bary splits each into 3 → 12
rr_bary_red = Adaptivity.compose_refinement_rules(rr_bary, rr_red)
@test Adaptivity.RefinementRuleType(rr_bary_red) isa Adaptivity.GenericRefinement
@test Adaptivity.get_polytope(rr_bary_red) == TRI
@test Adaptivity.num_subcells(rr_bary_red) == Adaptivity.num_subcells(rr_red) * Adaptivity.num_subcells(rr_bary)

# Red ∘ Bary: Bary splits into 3, Red splits each into 4 → 12
rr_red_bary = Adaptivity.compose_refinement_rules(rr_red, rr_bary)
@test Adaptivity.RefinementRuleType(rr_red_bary) isa Adaptivity.GenericRefinement
@test Adaptivity.get_polytope(rr_red_bary) == TRI
@test Adaptivity.num_subcells(rr_red_bary) == Adaptivity.num_subcells(rr_bary) * Adaptivity.num_subcells(rr_red)

############################################################################################
# compose_refinement_rules: varargs (Red ∘ Red ∘ Red = 4^3 = 64)
############################################################################################

rr_red3 = Adaptivity.compose_refinement_rules(rr_red, rr_red, rr_red)
@test Adaptivity.RefinementRuleType(rr_red3) isa Adaptivity.GenericRefinement
@test Adaptivity.num_subcells(rr_red3) == Adaptivity.num_subcells(rr_red)^3

############################################################################################
# compose_glues: two successive coarsenings of a simplexified voronoi mesh
############################################################################################

model1 = Geometry.voronoi(Geometry.simplexify(CartesianDiscreteModel((0,1,0,1),(3,3))))

pcells1 = Table([
  LinearIndices((4,4))[1:2,1:2],
  LinearIndices((4,4))[3:4,1:2],
  LinearIndices((4,4))[1:2,3:4],
  LinearIndices((4,4))[3:4,3:4],
])
ptopo1 = Geometry.PatchTopology(get_grid_topology(model1), pcells1)
model2, glue12 = Adaptivity.coarsen(model1, ptopo1; return_glue=true)

pcells2 = Table([[1,2],[3,4]])
ptopo2 = Geometry.PatchTopology(get_grid_topology(model2), pcells2)
model3, glue23 = Adaptivity.coarsen(model2, ptopo2; return_glue=true)

glue13 = Adaptivity.compose_glues(glue12, glue23)

Dc    = 2
topo1 = get_grid_topology(model1)
topo3 = get_grid_topology(model3)

# Face maps: length and zero-count checks for d = 0, 1, 2
for d in 0:Dc
  n2o = glue13.n2o_faces_map[d+1]
  @test length(n2o) == num_faces(topo1, d)
  if d == Dc
    @test count(iszero, n2o) == 0
  else
    @test count(!iszero, n2o) == num_faces(topo3, d)
  end
end

# Refinement rules: one per mesh3 cell, polytope matches glue23
n3 = num_cells(topo3)
@test length(glue13.refinement_rules) == n3
for rr in glue13.refinement_rules
  @test Adaptivity.get_polytope(rr) == Adaptivity.get_polytope(glue23.refinement_rules[1])
end

############################################################################################
# compress_adaptivity: 3-level uniform red refinement on TRI
############################################################################################

base_tri = Geometry.simplexify(CartesianDiscreteModel((0,1,0,1),(3,3)))
m4_tri   = refine(refine(refine(base_tri)))
c_tri    = Adaptivity.compress_adaptivity(m4_tri)

@test num_cells(c_tri.model) == num_cells(base_tri) * 4^3
@test !isa(c_tri.parent, AdaptedDiscreteModel)
@test length(c_tri.glue.refinement_rules) == num_cells(base_tri)
@test all(rr -> Adaptivity.num_subcells(rr) == 4^3, c_tri.glue.refinement_rules)

############################################################################################
# compress_adaptivity: 3-level uniform red refinement on QUAD
############################################################################################

base_quad = Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(3,3)))
m4_quad   = refine(refine(refine(base_quad)))
c_quad    = Adaptivity.compress_adaptivity(m4_quad)

@test num_cells(c_quad.model) == num_cells(base_quad) * 4^3
@test !isa(c_quad.parent, AdaptedDiscreteModel)
@test length(c_quad.glue.refinement_rules) == num_cells(base_quad)
@test all(rr -> Adaptivity.num_subcells(rr) == 4^3, c_quad.glue.refinement_rules)

end # module
