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
V_f = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
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
x_c = randn(num_free_dofs(U_c))
uh_c = FEFunction(U_c,x_c)
uh_f = change_domain(uh_c,trian)

# Coarse FEBasis -> Fine CellField
feb_c = get_fe_basis(V_c)
feb_f = get_fe_basis(V_f)
feb_c2f = change_domain(feb_c,trian)

# Coarse assembly
assem_c = SparseMatrixAssembler(U_c,V_c)
contr_c = bil(uh_c,feb_c,dΩ_c)
vecdata = collect_cell_vector(V_c,contr_c)
vec_c   = assemble_vector(assem_c,vecdata)

# Assembly of fine feb + c2f fefunc into Ω_f
assem_f = SparseMatrixAssembler(U_f,V_f)
contr_f = bil(uh_f,feb_f,dΩ_f)
vecdata = collect_cell_vector(V_f,contr_f)
vec_f   = assemble_vector(assem_f,vecdata)

# Assembly of c2f feb + c2f fefunc into Ω_f
assem_c2f = SparseMatrixAssembler(U_c,V_c)
contr_c2f = bil(uh_f,feb_c2f,dΩ_f)
contr_c2f2c = Refinement.merge_contr_cells(contr_c2f,trian,ctrian)
vecdata = collect_cell_vector(V_c,contr_c2f2c)
vec_c2f = assemble_vector(assem_c2f,vecdata)

@test vec_c ≈ vec_c2f

# Coarse FEFunction -> Fine FEFunction
op_c2f = RefinementTransferOperator(U_c,U_f)
x = copy(x_c)
y = zeros(num_free_dofs(U_f))
uh_c = FEFunction(U_c,x)
mul!(y,op_c2f,copy(x))
uh_f = FEFunction(U_f,y)

pts = map(x -> VectorValue(rand(2)),1:10)
v_c = map(p -> uh_c(p), pts)
v_f = map(p -> uh_f(p), pts)
@test v_c ≈ v_f

# Fine FEFunction -> Coarse FEFunction
op_f2c = RefinementTransferOperator(U_f,U_c)
x = randn(num_free_dofs(U_f))
y = zeros(num_free_dofs(U_c))
uh_f = FEFunction(U_f,x)
mul!(y,op_f2c,copy(x))
uh_c = FEFunction(U_c,y)

pts = map(x -> VectorValue(rand(2)),1:10)
v_c = map(p -> uh_c(p), pts)
v_f = map(p -> uh_f(p), pts)
@test v_c ≈ v_f

#end