module UniformRefinementTests
  
using Test
using Gridap
using Gridap.Adaptivity
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Adaptivity: uniformly_refine

using ..EdgeBasedRefinementTests: test_grid_transfers

# Setup base models
has_affine_map = true
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model1 = UnstructuredDiscreteModel(cart_model;has_affine_map)
model2 = simplexify(model1)

cart_model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
model3 = UnstructuredDiscreteModel(cart_model;has_affine_map)
model4 = simplexify(model3)

model5 = Geometry.DiscreteModelMock()

periodic_model = CartesianDiscreteModel((0,1,0,1),(4,4);isperiodic=(true,false))
model6 = UnstructuredDiscreteModel(periodic_model;has_affine_map)

n = 3

visualize = false
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
ref_model = refine(model4,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_tet_$n"))
test_grid_transfers(model4,ref_model,1)

# Mock
ref_model = refine(model5,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_mock_$n"))
test_grid_transfers(model5,ref_model,1)

# Periodic QUAD
ref_model = refine(model6,n)
visualize && writevtk(ref_model,joinpath(path,"uniform_periodic_quad_$n"))
# test_grid_transfers(model6,ref_model,1)

# Partial refinement
cell_refine_masks = [1,3,5]
ref_model = uniformly_refine(model5,n,cell_refine_masks)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_mock_partial_$n"))
test_grid_transfers(model5,ref_model,1)

end