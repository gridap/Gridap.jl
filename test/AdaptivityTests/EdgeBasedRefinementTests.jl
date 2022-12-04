module EdgeBasedRefinementTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model      = UnstructuredDiscreteModel(cart_model)

# Homogeneous refinement
ref_model1 = refine(model)
trian1 = Triangulation(ref_model1.model)

# Propagate to all-red
ref_model2 = refine(model;cells_to_refine=[1,6,11,16])
trian2 = Triangulation(ref_model2.model)

# Red-Green refinement
ref_model3 = refine(model;cells_to_refine=[1,6,16])
trian3 = Triangulation(ref_model3.model)

ref_model4 = refine(model;cells_to_refine=[6,7,10,11])
trian4 = Triangulation(ref_model4.model)

visualize = false
if visualize
  writevtk(trian1,"test/AdaptivityTests/red1_ref")
  writevtk(trian2,"test/AdaptivityTests/red2_ref")
  writevtk(trian3 ,"test/AdaptivityTests/redgreen1_ref")
  writevtk(trian4 ,"test/AdaptivityTests/redgreen2_ref")
end

end