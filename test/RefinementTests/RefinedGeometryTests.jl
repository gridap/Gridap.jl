module RefinedGeometryTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Refinement
using Gridap.ReferenceFEs
using FillArrays

function build_refined_model()
  D = (0,1,0,1)

  # Models
  parent = CartesianDiscreteModel(D,(1,1))
  child  = CartesianDiscreteModel(D,(2,2))

  # Glue 
  faces_map      = [Int[],Int[],Int[1,1,1,1]]
  fcell_child_id = Int[1,2,3,4]
  reffe          = LagrangianRefFE(Float64,QUAD,1)
  ref_cell_map   = get_f2c_ref_cell_map(reffe)
  glue = Gridap.Refinement.RefinementGlue(faces_map,fcell_child_id,ref_cell_map)

  # RefinedModel
  model = RefinedDiscreteModel(child,parent,glue)
  return model
end

# Get refined model and triangulation
model = build_refined_model()
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