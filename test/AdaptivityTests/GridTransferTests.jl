module GridTransferTests

using Test
using Gridap
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

# Solutions and weakform
sol(x) = x[1] + x[2]
bil(uh,vh,dΩ) = ∫(uh⋅vh)*dΩ
solver = BackslashSolver()

# Get refined model and triangulation
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model = refine(cart_model; num_refinements=2)
trian = Triangulation(model)

# Triangulations
ftrian = Triangulation(get_model(model))
ctrian = Triangulation(get_parent(model))

# Measures
dΩ_f   = Measure(ftrian,2)
dΩ_c   = Measure(ctrian,2)
dΩ_cf  = Measure(ctrian,trian,1)

cell_quad = Gridap.CellData.get_cell_quadrature(dΩ_cf)
dΩ_cf_bis = Measure(ctrian,trian,cell_quad)

# FESpaces
reffe = ReferenceFE(lagrangian,Float64,1)
V_f = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
U_f = TrialFESpace(V_f,sol)
V_c = TestFESpace(get_parent(model),reffe;conformity=:H1,dirichlet_tags="boundary")
U_c = TrialFESpace(V_c,sol)

# CellField: Coarse -> Fine
cf_c_phy = CellField(sol,ctrian)
cf_c_ref = change_domain(cf_c_phy,PhysicalDomain(),ReferenceDomain())
cf_f_ref_ref = change_domain(cf_c_ref, trian, ReferenceDomain())
cf_f_ref_phy = change_domain(cf_c_ref, trian, PhysicalDomain())
cf_f_phy_ref = change_domain(cf_c_phy, trian, ReferenceDomain())
cf_f_phy_phy = change_domain(cf_c_phy, trian, PhysicalDomain())

pts = map(x -> VectorValue(rand(2)),1:10)
v_r = map(p -> sol(p) , pts)
v_c = map(p -> cf_c_ref(p), pts)
v_f_ref_ref = map(p -> cf_f_ref_ref(p), pts)
v_f_ref_phy = map(p -> cf_f_ref_phy(p), pts)
v_f_phy_ref = map(p -> cf_f_phy_ref(p), pts)
v_f_phy_phy = map(p -> cf_f_phy_phy(p), pts)
@test v_r ≈ v_c
@test v_r ≈ v_f_ref_ref
@test v_r ≈ v_f_ref_phy
@test v_r ≈ v_f_phy_ref
@test v_r ≈ v_f_phy_phy

# CellField: Fine -> Coarse
cf_f_phy = CellField(sol,trian)
cf_f_ref = change_domain(cf_f_phy,PhysicalDomain(),ReferenceDomain())
cf_c_ref_ref = change_domain(cf_f_ref, ctrian, ReferenceDomain())
cf_c_ref_phy = change_domain(cf_f_ref, ctrian, PhysicalDomain())
cf_c_phy_ref = change_domain(cf_f_phy, ctrian, ReferenceDomain())
cf_c_phy_phy = change_domain(cf_f_phy, ctrian, PhysicalDomain())

pts = map(x -> VectorValue(rand(2)),1:10)
v_r = map(p -> sol(p) , pts)
v_f = map(p -> cf_f_ref(p), pts)
v_c_ref_ref = map(p -> cf_c_ref_ref(p), pts)
v_c_ref_phy = map(p -> cf_c_ref_phy(p), pts)
v_c_phy_ref = map(p -> cf_c_phy_ref(p), pts)
v_c_phy_phy = map(p -> cf_c_phy_phy(p), pts)
@test v_r ≈ v_f
@test v_r ≈ v_c_ref_ref
@test v_r ≈ v_c_ref_phy
@test v_r ≈ v_c_phy_ref
@test v_r ≈ v_c_phy_phy

# Coarse FEFunction -> Fine CellField
uh_c = interpolate(sol,U_c)
cf_f2 = change_domain(uh_c,trian,ReferenceDomain())
v_f2 = map(p -> cf_f2(p), pts)
@test v_r ≈ v_f2

# Coarse FEBasis -> Fine CellField
feb_c = get_fe_basis(V_c)
feb_c2f = change_domain(feb_c,trian,ReferenceDomain())

# Coarse FEFunction -> Fine FEFunction, by interpolation
uh_f_inter = interpolate(uh_c,U_f)

v_f_inter = map(p -> uh_f_inter(p), pts)
@test v_r ≈ v_f_inter

# Fine FEFunction -> Coarse FEFunction, by interpolation
uh_f = interpolate(sol,U_f)
uh_c_inter = interpolate(uh_f,U_c)

v_c_inter = map(p -> uh_c_inter(p), pts)
@test v_r ≈ v_c_inter

# Coarse FEFunction -> Fine FEFunction, by projection
af(u,v)  = ∫(v⋅u)*dΩ_f
lf(v)    = ∫(v⋅uh_c)*dΩ_f
opf      = AffineFEOperator(af,lf,U_f,V_f)
uh_f_pr = solve(opf)

v_f_pr = map(p -> uh_f_pr(p), pts)
@test v_r ≈ v_f_pr

# Fine FEFunction -> Coarse FEFunction, by projection
ac(u,v) = ∫(v⋅u)*dΩ_c
lc(v)   = ∫(v⋅uh_f_inter)*dΩ_cf
opc     = AffineFEOperator(ac,lc,U_c,V_c)
uh_c_pr = solve(opc)

v_c_pr = map(p -> uh_c_pr(p), pts)
@test v_c_pr ≈ v_r

end