module EdgeBasedRefinementTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

visualize = false

# Refining meshes of QUADs
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model      = UnstructuredDiscreteModel(cart_model)

## Homogeneous refinement
ref_model1 = refine(model)
trian1 = Triangulation(ref_model1.model)
visualize && writevtk(trian1,"test/AdaptivityTests/ref_model1")

## Propagate to all-red
ref_model2 = refine(model;cells_to_refine=[1,6,11,16])
trian2 = Triangulation(ref_model2.model)
visualize && writevtk(trian2,"test/AdaptivityTests/ref_model2")

## Red-Green refinement
ref_model3 = refine(model;cells_to_refine=[1,6,16])
trian3 = Triangulation(ref_model3.model)
visualize && writevtk(trian3,"test/AdaptivityTests/ref_model3")

ref_model4 = refine(model;cells_to_refine=[6,7,10,11])
trian4 = Triangulation(ref_model4.model)
visualize && writevtk(trian4,"test/AdaptivityTests/ref_model4")

# Refining meshes of TRIans
model2 = simplexify(model)
visualize && writevtk(Triangulation(model2),"test/AdaptivityTests/base_model2")

ref_model5 = refine(model2)
trian5 = Triangulation(ref_model5.model)
visualize && writevtk(trian5,"test/AdaptivityTests/ref_model5")

ref_model6 = refine(model2;cells_to_refine=[1,6,16])
trian6 = Triangulation(ref_model6.model)
visualize && writevtk(trian6,"test/AdaptivityTests/ref_model6")

# Testing FESpaces
sol(x) = x[1] + x[2]
reffe = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(ref_model6,reffe;conformity=:H1,dirichlet_tags="boundary")
U = TrialFESpace(V,sol)

end