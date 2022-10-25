module RefinedGeometryTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Refinement
using Gridap.ReferenceFEs
using FillArrays

# Get refined model and triangulation
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model  = refine(cart_model; num_refinements=2)
test_discrete_model(model)
trian = Triangulation(model)
test_triangulation(trian)
test_triangulation(trian.trian)
@test isa(trian, RefinedTriangulation)

# Get members
fmodel = get_model(model)
cmodel = get_parent(model)
glue   = get_refinement_glue(model)

# Triangulations
ftrian = Triangulation(fmodel)
ctrian = Triangulation(cmodel)

# Choosing targets
model2 = refine(model; num_refinements=2)
trian2 = Triangulation(model2)
@test best_target(trian,ftrian) === trian
@test best_target(trian,ctrian) === trian
@test best_target(trian,trian2) === trian2

end