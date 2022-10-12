#module RefinedDiscreteModelsTests

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

# CellField: Coarse -> Fine
func(x) = x[1] + x[2]
cf_c = change_domain(CellField(func,ctrian),PhysicalDomain(),ReferenceDomain())
cf_f = change_domain_c2f(cf_c, ftrian, model)

pts = map(x->VectorValue(rand(2)),1:10)
v_r = map(p-> func(p), pts)  # Real values
v_c = map(p -> cf_c(p), pts) # Values by Coarse CellField
v_f = map(p -> cf_f(p), pts) # Values by Fine CellField
@test v_r ≈ v_c
@test v_r ≈ v_f

#end