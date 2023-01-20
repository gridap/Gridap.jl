module EdgeBasedRefinementTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

"""
  If OrientationStyle(topo)==true, checks if the topology is indeed oriented.
  If OrientationStyle(topo)==false or the topology is not oriented, returns false.
"""
function test_grid_topology_orientation(topo::GridTopology{Dc,Dp}) where {Dc,Dp}
  orientation = OrientationStyle(topo)
  isa(orientation,NonOriented) && (return false)

  nC = num_faces(topo,Dc)
  c2n_map     = get_faces(topo,Dc,0)
  is_oriented = true; iC = 1
  while(is_oriented && iC <= nC)
    nodes = c2n_map[iC]
    is_oriented = (is_oriented && issorted(nodes))
    iC += 1
  end
  return is_oriented
end

visualize = false

# Refining meshes of QUADs
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model1     = UnstructuredDiscreteModel(cart_model)

## Homogeneous refinement
ref_model1 = refine(model1)
trian1 = Triangulation(ref_model1.model)
visualize && writevtk(trian1,"test/AdaptivityTests/ref_model1")

## Propagate to all-red
ref_model2 = refine(model1;cells_to_refine=[1,6,11,16])
trian2 = Triangulation(ref_model2.model)
visualize && writevtk(trian2,"test/AdaptivityTests/ref_model2")

## Red-Green refinement
ref_model3 = refine(model1;cells_to_refine=[1,6,16])
trian3 = Triangulation(ref_model3.model)
visualize && writevtk(trian3,"test/AdaptivityTests/ref_model3")

ref_model4 = refine(model1;cells_to_refine=[6,7,10,11])
trian4 = Triangulation(ref_model4.model)
visualize && writevtk(trian4,"test/AdaptivityTests/ref_model4")

# Refining meshes of TRIans
model2 = simplexify(model1)
visualize && writevtk(Triangulation(model2),"test/AdaptivityTests/base_model2")

ref_model5 = refine(model2)
trian5 = Triangulation(ref_model5.model)
visualize && writevtk(trian5,"test/AdaptivityTests/ref_model5")

ref_model6 = refine(model2;cells_to_refine=[1,6,16])
trian6 = Triangulation(ref_model6.model)
visualize && writevtk(trian6,"test/AdaptivityTests/ref_model6")

# Testing FE related functionality
parent = model2
model  = ref_model6

topo = get_grid_topology(model)
@test isa(OrientationStyle(topo),Oriented)
@test test_grid_topology_orientation(topo) == true

u((x,y)) = 2*VectorValue(-y,x)
reffe = ReferenceFE(nedelec,Float64,1)
Vh = TestFESpace(model,reffe;dirichlet_tags="boundary")
Uh = TrialFESpace(Vh,u)
VH = TestFESpace(parent,reffe;dirichlet_tags="boundary")
UH = TrialFESpace(VH,u)

uh = interpolate(u,Uh)
uH = interpolate(u,UH)

Ωh = Triangulation(model)
dΩh = Measure(Ωh,3)

e = u - uH
el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
@test el2 < 1.0e-10

end