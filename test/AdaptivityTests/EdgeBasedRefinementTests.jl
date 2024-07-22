module EdgeBasedRefinementTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

function l2_error(u1,u2,dΩ)
  eh = u1-u2
  return sum(∫(eh⋅eh)*dΩ)
end

function test_grid_transfers(D,parent,model,order)
  sol(x) = sum(x)
  qorder = 2*order+1

  # Triangulations
  trian = Triangulation(model)
  ctrian = Triangulation(parent)

  glue = get_adaptivity_glue(model)
  rrules = Adaptivity.get_old_cell_refinement_rules(glue)

  # Measures
  dΩ_f  = Measure(trian,qorder)
  dΩ_c  = Measure(ctrian,qorder)
  dΩ_cf = Measure(ctrian,trian,qorder)

  cell_quad = Gridap.CellData.get_cell_quadrature(dΩ_cf)
  dΩ_cf_bis = Measure(ctrian,trian,cell_quad)

  # FESpaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_f = TestFESpace(model,reffe;dirichlet_tags="boundary")
  U_f = TrialFESpace(V_f,sol)
  V_c = TestFESpace(parent,reffe;dirichlet_tags="boundary")
  U_c = TrialFESpace(V_c,sol)
  V_c_fast = TestFESpace(parent,rrules,reffe;dirichlet_tags="boundary")
  U_c_fast = TrialFESpace(V_c_fast,sol)

  # CellField: Coarse -> Fine
  cf_c_phy = CellField(sol,ctrian)
  cf_c_ref = change_domain(cf_c_phy,PhysicalDomain(),ReferenceDomain())
  cf_f_ref_ref = change_domain(cf_c_ref, trian, ReferenceDomain())
  cf_f_ref_phy = change_domain(cf_c_ref, trian, PhysicalDomain())
  cf_f_phy_ref = change_domain(cf_c_phy, trian, ReferenceDomain())
  cf_f_phy_phy = change_domain(cf_c_phy, trian, PhysicalDomain())

  # CellField: Fine -> Coarse
  cf_f_phy = CellField(sol,trian)
  cf_f_ref = change_domain(cf_f_phy,PhysicalDomain(),ReferenceDomain())
  cf_c_ref_ref = change_domain(cf_f_ref, ctrian, ReferenceDomain())
  cf_c_ref_phy = change_domain(cf_f_ref, ctrian, PhysicalDomain())
  cf_c_phy_ref = change_domain(cf_f_phy, ctrian, ReferenceDomain())
  cf_c_phy_phy = change_domain(cf_f_phy, ctrian, PhysicalDomain())

  @test l2_error(cf_f_ref_ref,cf_f_ref,dΩ_f) < 1.e-8
  @test l2_error(cf_f_ref_phy,cf_f_phy,dΩ_f) < 1.e-8
  @test l2_error(cf_f_phy_ref,cf_f_ref,dΩ_f) < 1.e-8
  @test l2_error(cf_f_phy_phy,cf_f_phy,dΩ_f) < 1.e-8

  @test l2_error(cf_c_ref_ref,cf_c_ref,dΩ_c) < 1.e-8
  @test l2_error(cf_c_ref_phy,cf_c_phy,dΩ_c) < 1.e-8
  @test l2_error(cf_c_phy_ref,cf_c_ref,dΩ_c) < 1.e-8
  @test l2_error(cf_c_phy_phy,cf_c_phy,dΩ_c) < 1.e-8

  # Coarse FEFunction -> Fine FEFunction, by interpolation
  uh_c = interpolate(sol,U_c)
  uh_f_inter  = interpolate(uh_c,U_f)
  uh_f_inter2 = interpolate_everywhere(uh_c,U_f)
  uh_f_inter3 = interpolate_dirichlet(uh_c,U_f)

  # Fine FEFunction -> Coarse FEFunction, by interpolation
  uh_f = interpolate(sol,U_f)
  uh_c_inter  = interpolate(uh_f,U_c)
  uh_c_inter2 = interpolate_everywhere(uh_f,U_c)
  uh_c_inter3 = interpolate_dirichlet(uh_f,U_c)
  uh_c_inter4  = interpolate(uh_f,U_c_fast)

  @test l2_error(uh_f,uh_f_inter,dΩ_f) < 1.e-8
  @test l2_error(uh_f,uh_f_inter2,dΩ_f) < 1.e-8

  @test l2_error(uh_c,uh_c_inter,dΩ_c) < 1.e-8
  @test l2_error(uh_c,uh_c_inter2,dΩ_c) < 1.e-8
  @test l2_error(uh_c,uh_c_inter4,dΩ_c) < 1.e-8

  # Coarse FEFunction -> Fine FEFunction, by projection
  af(u,v)  = ∫(v⋅u)*dΩ_f
  lf(v)    = ∫(v⋅uh_c)*dΩ_f
  opf      = AffineFEOperator(af,lf,U_f,V_f)
  uh_f_pr = solve(opf)
  @test l2_error(uh_f,uh_f_pr,dΩ_f) < 1.e-8

  # Fine FEFunction -> Coarse FEFunction, by projection
  ac(u,v) = ∫(v⋅u)*dΩ_c
  lc(v)   = ∫(v⋅uh_f_inter)*dΩ_cf
  opc     = AffineFEOperator(ac,lc,U_c,V_c)
  uh_c_pr = solve(opc)
  @test l2_error(uh_c,uh_c_pr,dΩ_c) < 1.e-8
end

visualize = false

############################################################################################
### Red-Green refinement

## A) 2D meshes - QUADs

cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model1     = UnstructuredDiscreteModel(cart_model)

# Homogeneous refinement
ref_model1 = refine(model1)
trian1 = Triangulation(ref_model1.model)
visualize && writevtk(trian1,"test/AdaptivityTests/ref_model1")
test_grid_transfers(2,model1,ref_model1,1)

# Propagate to all-red
ref_model2 = refine(model1;cells_to_refine=[1,6,11,16])
trian2 = Triangulation(ref_model2.model)
visualize && writevtk(trian2,"test/AdaptivityTests/ref_model2")
test_grid_transfers(2,model1,ref_model2,1)

# Red-Green refinement
ref_model3 = refine(model1;cells_to_refine=[1,6,16])
trian3 = Triangulation(ref_model3.model)
visualize && writevtk(trian3,"test/AdaptivityTests/ref_model3")
#test_grid_transfers(2,model1,ref_model3,1)

ref_model4 = refine(model1;cells_to_refine=[6,7,10,11])
trian4 = Triangulation(ref_model4.model)
visualize && writevtk(trian4,"test/AdaptivityTests/ref_model4")
#test_grid_transfers(2,model1,ref_model4,1)

## B) 2D meshes - TRIs

model2 = simplexify(model1)
visualize && writevtk(Triangulation(model2),"test/AdaptivityTests/base_model2")

ref_model5 = refine(model2)
trian5 = Triangulation(ref_model5.model)
visualize && writevtk(trian5,"test/AdaptivityTests/ref_model5")
test_grid_transfers(2,model2,ref_model5,1)

ref_model6 = refine(model2;cells_to_refine=[1,6,16])
trian6 = Triangulation(ref_model6.model)
visualize && writevtk(trian6,"test/AdaptivityTests/ref_model6")
test_grid_transfers(2,model2,ref_model6,1)

## C) 3D meshes - HEXs

cart_model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
model3     = UnstructuredDiscreteModel(cart_model)

ref_model7 = refine(model3)
trian7 = Triangulation(ref_model7.model)
visualize && writevtk(trian7,"test/AdaptivityTests/ref_model7")
test_grid_transfers(3,model3,ref_model7,1)

## D) 3D meshes - TETs

model4 = simplexify(model3)
ref_model8 = refine(model4)
trian8 = Triangulation(ref_model8.model)
visualize && writevtk(trian8,"test/AdaptivityTests/ref_model8")
test_grid_transfers(3,model4,ref_model8,1)

############################################################################################
### Newest Vertex Bisection refinement (longest edge bisection)

# Refine all edges using NVB
ref_model9 = refine(model2, refinement_method = "nvb")
trian9 = Triangulation(ref_model9.model)
visualize && writevtk(trian9, "test/AdaptivityTests/ref_model7")
test_grid_transfers(2, model2, ref_model9, 1)

# Refine all edges using NVB
# Mark edges such that blue, and double_blue refinement are triggered
ref_model10 = refine(model2, refinement_method = "nvb", cells_to_refine = [4, 9])
trian10 = Triangulation(ref_model10)
visualize && writevtk(trian10, "test/AdaptivityTests/ref_model8")
ref_model11 = refine(ref_model10, refinement_method = "nvb", cells_to_refine = [1, 3, 4, 11])
trian11 = Triangulation(ref_model11.model)
visualize && writevtk(trian11, "test/AdaptivityTests/ref_model9")
test_grid_transfers(2, ref_model10, ref_model11, 1)

############################################################################################
### Barycentric refinement

## A) 2D meshes - TRIs

ref_model12 = refine(model2, refinement_method = "barycentric")
trian12 = Triangulation(ref_model12.model)
visualize && writevtk(trian12, "test/AdaptivityTests/ref_model12")
test_grid_transfers(2, model2, ref_model12, 1)

ref_model13 = refine(model2, refinement_method = "barycentric", cells_to_refine = [1, 6, 8])
trian13 = Triangulation(ref_model13.model)
visualize && writevtk(trian13, "test/AdaptivityTests/ref_model13")
test_grid_transfers(2, model2, ref_model13, 1)

## B) 3D meshes - TETs

ref_model14 = refine(model4, refinement_method = "barycentric")
trian14 = Triangulation(ref_model14.model)
visualize && writevtk(trian13, "test/AdaptivityTests/ref_model13")
test_grid_transfers(3, model4, ref_model14, 1)

ref_model15 = refine(model4, refinement_method = "barycentric", cells_to_refine = [1, 6, 8])
trian15 = Triangulation(ref_model15.model)
visualize && writevtk(trian15, "test/AdaptivityTests/ref_model15")
test_grid_transfers(3, model4, ref_model15, 1)

############################################################################################
### Simplexify refinement

## A) 2D meshes - QUADs
ref_model16 = refine(model1, refinement_method = "simplexify")
trian16 = Triangulation(ref_model16.model)
visualize && writevtk(trian16, "test/AdaptivityTests/ref_model16")
test_grid_transfers(2, model1, ref_model16, 1)

## B) 3D meshes - HEXs
ref_model17 = refine(model3, refinement_method = "simplexify")
trian17 = Triangulation(ref_model17.model)
visualize && writevtk(trian17, "test/AdaptivityTests/ref_model17")
test_grid_transfers(3, model3, ref_model17, 1)

## C) 2D meshes - QUADs (Positive Volume)
ref_model18 = refine(model1, refinement_method = "simplexify", positive = true)
trian18 = Triangulation(ref_model18.model)
visualize && writevtk(trian18, "test/AdaptivityTests/ref_model18")
test_grid_transfers(2, model1, ref_model18, 1)

## D) 3D meshes - HEXs (Positive Volume)
ref_model19 = refine(model3, refinement_method = "simplexify", positive = true)
trian19 = Triangulation(ref_model19.model)
visualize && writevtk(trian19, "test/AdaptivityTests/ref_model19")
test_grid_transfers(3, model3, ref_model19, 1)
end
