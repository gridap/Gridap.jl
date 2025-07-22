module UnstructuredUniformRefinementTests
  
using Test
using Gridap
using Gridap.Adaptivity
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Adaptivity: test_unstructured_uniform_refinement

using .EdgeBasedRefinementTests: test_grid_transfers

test_unstructured_uniform_refinement()

# Setup base models

cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model1 = UnstructuredDiscreteModel(cart_model)
model2 = simplexify(model1)

cart_model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
model3 = UnstructuredDiscreteModel(cart_model)
model4 = simplexify(model3)

model5 = Geometry.DiscreteModelMock()

n = 3

visualize = true
if visualize
  path = mkpath("tmp/")
end

# QUAD
ref_model = refine(model1,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_quad_$n"))
test_grid_transfers(model1,ref_model,1)

# TRI
ref_model = refine(model2,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_tri_$n"))
test_grid_transfers(model2,ref_model,1)

# HEX
ref_model = refine(model3,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_hex_$n"))
test_grid_transfers(model3,ref_model,1)

# TET
# bug
ref_model = refine(model4,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_tet_$n"))
test_grid_transfers(model4,ref_model,1)

# Mock
# bug in calling TestFESpace(model5,rrules,ReferenceFE(lagrangian,Float64,1))
ref_model = refine(model5,n)
# visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_mock_$n"))
glue = get_adaptivity_glue(ref_model)
rrules = Adaptivity.get_old_cell_refinement_rules(glue)
TestFESpace(model5,rrules,ReferenceFE(lagrangian,Float64,1))
# test_grid_transfers(model5,ref_model,1)

end