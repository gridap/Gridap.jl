module RefinedGeometryTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Refinement
using Gridap.ReferenceFEs
using FillArrays

# Get refined model and triangulation
model = RefinedCartesianDiscreteModel((0,1,0,1),4,2)
test_discrete_model(model)
trian = Triangulation(model)
test_triangulation(trian)
test_triangulation(trian.trian)
@test isa(trian, RefinedTriangulation)

# Get members
fmodel = get_model(model)
cmodel = get_parent(model)
glue   = get_glue(model)

# Triangulations
ftrian = Triangulation(fmodel)
ctrian = Triangulation(cmodel)

end