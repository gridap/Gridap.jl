module EdgeBasedRefinementTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

cart_model = CartesianDiscreteModel((0,1,0,1),(2,2))
model      = UnstructuredDiscreteModel(cart_model)


ref_model1 = refine(model)
trian1 = Triangulation(ref_model1.model)
writevtk(trian1,"red1_ref")

ref_model2 = refine(model;cells_to_refine=[1,4])
trian2 = Triangulation(ref_model2.model)
writevtk(trian2,"red2_ref")

ref_model3 = refine(model;cells_to_refine=[1])
trian3 = Triangulation(ref_model3.model)
writevtk(trian3 ,"redgreen1_ref")

ref_model4 = refine(model;cells_to_refine=[2,4])
trian4 = Triangulation(ref_model4.model)
writevtk(trian3 ,"redgreen2_ref")