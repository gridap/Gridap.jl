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

function test_grid_transfers(parent,model,order)
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

# Setup base models

cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model1 = UnstructuredDiscreteModel(cart_model)
model2 = simplexify(model1)

cart_model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
model3 = UnstructuredDiscreteModel(cart_model)
model4 = simplexify(model3)

visualize = false
if visualize
  path = mkpath("tmp/")
end

############################################################################################
### Red-Green refinement

## A) 2D meshes - QUADs

# Homogeneous refinement
ref_model = refine(model1)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen1"))
test_grid_transfers(model1,ref_model,1)

# Propagate to all-red
ref_model = refine(model1;cells_to_refine=[1,6,11,16])
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen2"))
test_grid_transfers(model1,ref_model,1)

# Red-Green refinement
ref_model = refine(model1;cells_to_refine=[1,6,16])
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen3"))
#test_grid_transfers(model1,ref_model,1)

ref_model = refine(model1;cells_to_refine=[6,7,10,11])
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen4"))
#test_grid_transfers(model1,ref_model,1)

## B) 2D meshes - TRIs

ref_model = refine(model2)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen5"))
test_grid_transfers(model2,ref_model,1)

ref_model = refine(model2;cells_to_refine=[1,6,16])
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen6"))
test_grid_transfers(model2,ref_model,1)

## C) 3D meshes - HEXs

ref_model = refine(model3)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen7"))
test_grid_transfers(model3,ref_model,1)

## D) 3D meshes - TETs

ref_model = refine(model4)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"redgreen8"))
test_grid_transfers(model4,ref_model,1)

############################################################################################
### Newest Vertex Bisection refinement (longest edge bisection)

# Refine all edges using NVB
ref_model = refine(model2, refinement_method = "nvb")
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"nvb1"))
test_grid_transfers(model2, ref_model, 1)

# Refine all edges using NVB
# Mark edges such that blue, and double_blue refinement are triggered
_ref_model = refine(model2, refinement_method = "nvb", cells_to_refine = [4, 9])
visualize && writevtk(Triangulation(_ref_model.model),joinpath(path,"nvb2"))
ref_model = refine(_ref_model, refinement_method = "nvb", cells_to_refine = [1, 3, 4, 11])
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"nvb3"))
test_grid_transfers(_ref_model, ref_model, 1)

############################################################################################
### Barycentric refinement

## A) 2D meshes - TRIs

ref_model = refine(model2, refinement_method = "barycentric")
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"barycentric1"))
test_grid_transfers(model2, ref_model, 1)

ref_model = refine(model2, refinement_method = "barycentric", cells_to_refine = [1, 6, 8])
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"barycentric2"))
test_grid_transfers(model2, ref_model, 1)

## B) 3D meshes - TETs

ref_model = refine(model4, refinement_method = "barycentric")
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"barycentric3"))
test_grid_transfers(model4, ref_model, 1)

ref_model = refine(model4, refinement_method = "barycentric", cells_to_refine = [1, 6, 8])
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"barycentric4"))
test_grid_transfers(model4, ref_model, 1)

## C) 2D meshes - QUADs

ref_model = refine(model1, refinement_method = "barycentric")
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"barycentric5"))
test_grid_transfers(model1, ref_model, 1)

## D) 3D meshes - HEXs

ref_model = refine(model3, refinement_method = "barycentric")
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"barycentric6"))
test_grid_transfers(model3, ref_model, 1)

############################################################################################
### Simplexify refinement

## A) 2D meshes - QUADs
ref_model = refine(model1, refinement_method = "simplexify")
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"simplexify1"))
test_grid_transfers(model1, ref_model, 1)

## B) 3D meshes - HEXs
ref_model = refine(model3, refinement_method = "simplexify")
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"simplexify2"))
test_grid_transfers(model3, ref_model, 1)

## C) 2D meshes - QUADs (Positive Volume)
ref_model = refine(model1, refinement_method = "simplexify", positive = true)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"simplexify3"))
test_grid_transfers(model1, ref_model, 1)

## D) 3D meshes - HEXs (Positive Volume)
ref_model = refine(model3, refinement_method = "simplexify", positive = true)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"simplexify4"))
test_grid_transfers(model3, ref_model, 1)

end
