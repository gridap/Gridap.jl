module CartesianRefinementTests

using Test
using Gridap
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

function l2_error(u1,u2,dΩ)
  eh = u1-u2
  return sum(∫(eh⋅eh)*dΩ)
end

function test_grid_transfers(D,model,parent,order)
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

for D = 1:3
  order  = 1
  domain = Tuple(repeat([0,1],D))

  cart_model = CartesianDiscreteModel(domain,Tuple(fill(4,D)))
  model1 = refine(cart_model,2)
  model2 = refine(model1,2)
  
  test_grid_transfers(D,model1,cart_model,order)
  test_grid_transfers(D,model2,model1,order)
end

end