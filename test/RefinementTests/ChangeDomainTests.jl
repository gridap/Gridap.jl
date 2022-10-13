#module ChangeDomainTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Refinement
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

sol(x) = x[1] + x[2]
bil(uh,vh,dΩ) = ∫(uh⋅vh)*dΩ

# Get refined model and triangulation
model = RefinedCartesianDiscreteModel((0,1,0,1),4,2)
trian = Triangulation(model)

# Triangulations
ftrian = Triangulation(get_model(model))
ctrian = Triangulation(get_parent(model))
dΩ_f   = Measure(ftrian,2)
dΩ_c   = Measure(ctrian,2)

# FESpaces
reffe = ReferenceFE(lagrangian,Float64,1)
V_f = TestFESpace(get_model(model),reffe;conformity=:H1,dirichlet_tags="boundary")
U_f = TrialFESpace(V_f,sol)
V_c = TestFESpace(get_parent(model),reffe;conformity=:H1,dirichlet_tags="boundary")
U_c = TrialFESpace(V_c,sol)

# CellField: Coarse -> Fine
cf_c = change_domain(CellField(sol,ctrian),PhysicalDomain(),ReferenceDomain())
cf_f = change_domain(cf_c, trian)

pts = map(x -> VectorValue(rand(2)),1:10)
v_r = map(p -> sol(p) , pts) # Real values
v_c = map(p -> cf_c(p), pts) # Values by Coarse CellField
v_f = map(p -> cf_f(p), pts) # Values by Fine CellField
@test v_r ≈ v_c
@test v_r ≈ v_f

# Coarse FEFunction -> Fine CellField
uh_c = FEFunction(U_c,randn(num_free_dofs(U_c)))
uh_f = change_domain(uh_c,trian)

# Coarse FEBasis -> Fine CellField
feb_c = get_fe_basis(V_c)
feb_f = change_domain(feb_c,trian)

# Check assembly
assem   = SparseMatrixAssembler(U_c,V_c)
contr_c = bil(uh_c,feb_c,dΩ_c)
contr_f = bil(uh_f,feb_f,dΩ_f)

vecdata = collect_cell_vector(V_c,contr_c)
vec_c = assemble_vector(assem,vecdata)

vecdata = collect_cell_vector(V_c,contr_f)
vec_f = assemble_vector(assem,vecdata)

display(vec_c)
display(vec_f)

# CoarseFEFunction -> Fine FEFunction


#end