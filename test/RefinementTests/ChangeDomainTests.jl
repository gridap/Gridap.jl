#module ChangeDomainTests

using Test
using Gridap
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Refinement
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

# Solutions and weakform
sol(x) = x[1] + x[2]
bil(uh,vh,dΩ) = ∫(uh⋅vh)*dΩ

"""
# Conjugate Gradient wrapper
struct CGLinearSolver <: LinearSolver end

struct CGSymbolicSetup{A} <: SymbolicSetup
  solver::A
end

struct CGNumericalSetup{A,B,C} <: NumericalSetup 
  ss :: A
  mat:: B
  caches:: C
end

function Gridap.Algebra.symbolic_setup(s::CGLinearSolver,mat::AbstractMatrix)
  CGSymbolicSetup(s)
end

function Gridap.Algebra.numerical_setup(ss::CGSymbolicSetup,mat::AbstractMatrix)
  caches = CG_get_caches(mat)
  CGNumericalSetup(ss,mat,caches)
end

function CG_get_caches(mat)
 return nothing
end

function Gridap.Algebra.solve!(x::AbstractVector,ns::CGNumericalSetup,b::AbstractVector)
  cg!(x,ns.mat,b)
end
"""

solver = BackslashSolver()

# Get refined model and triangulation
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model = refine(cart_model; num_refinements=2)
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
cf_f = change_domain(cf_c, trian, ReferenceDomain())

pts = map(x -> VectorValue(rand(2)),1:10)
v_r = map(p -> sol(p) , pts) # Real values
v_c = map(p -> cf_c(p), pts) # Values by Coarse CellField
v_f = map(p -> cf_f(p), pts) # Values by Fine CellField
@test v_r ≈ v_c
@test v_r ≈ v_f

# Coarse FEFunction -> Fine CellField
x_c = randn(num_free_dofs(U_c))
uh_c = FEFunction(U_c,x_c)
uh_f = change_domain(uh_c,trian,ReferenceDomain())

# Coarse FEBasis -> Fine CellField
feb_c = get_fe_basis(V_c)
feb_f = get_fe_basis(V_f)
feb_c2f = change_domain(feb_c,trian,ReferenceDomain())

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

# Assembly of c2f feb + fine fefunc into Ω_c
assem_c2f = SparseMatrixAssembler(U_c,V_c)
contr_c2f = bil(uh_f,feb_c2f,dΩ_f)
contr_c2f2c = Refinement.merge_contr_cells(contr_c2f,trian,ctrian)
vecdata = collect_cell_vector(V_c,contr_c2f2c)
vec_c2f = assemble_vector(assem_c2f,vecdata)

@test vec_c ≈ vec_c2f

# Coarse FEFunction -> Fine FEFunction using ProjectionTransferOperator
op_c2f = ProjectionTransferOperator(U_c,U_f;solver=solver)
y_f = zeros(num_free_dofs(U_f))
uh_c = FEFunction(U_c,x_c)
mul!(y_f,op_c2f,copy(x_c))
uh_f = FEFunction(U_f,y_f)

pts = map(x -> VectorValue(rand(2)),1:10)
v_c = map(p -> uh_c(p), pts)
v_f = map(p -> uh_f(p), pts)
@test v_c ≈ v_f

# Fine FEFunction -> Coarse FEFunction using ProjectionTransferOperator
op_f2c = ProjectionTransferOperator(U_f,U_c;solver=solver)
x_f = copy(y_f)
y_c = zeros(num_free_dofs(U_c))
uh_f = FEFunction(U_f,x_f)
mul!(y_c,op_f2c,copy(x_f))
uh_c = FEFunction(U_c,y_c)
@test y_c ≈ x_c

pts = map(x -> VectorValue(rand(2)),1:10)
v_c = map(p -> uh_c(p), pts)
v_f = map(p -> uh_f(p), pts)
@test v_c ≈ v_f

# Coarse to Fine, using interpolate()
uh_f_inter = interpolate(uh_c,U_f)

# Fine to Coarse, using interpolate()
uh_c_inter = interpolate(uh_f,U_c)

# Integrating with high-level API Coarse to Fine
dΩ_cf = CompositeMeasure(ctrian,trian,1)
vh_c = get_fe_basis(U_c)
vh_f = get_fe_basis(U_f)

integrand_fc = vh_f ⋅ uh_c
contr_fc_f = ∫(integrand_fc)*dΩ_f
vecdata = collect_cell_vector(U_f,contr_fc_f)
vec_f = assemble_vector(assem_f,vecdata)

integrand_cf = vh_c ⋅ uh_f
contr_cf_f = ∫(integrand_cf)*dΩ_f
vecdata = collect_cell_vector(U_f,contr_cf_f)

# Integrating with high-level API Fine to Coarse
contr_fc_c = ∫(integrand_fc)*dΩ_cf
vecdata = collect_cell_vector(U_c,contr_fc_c)

contr_cf_c = ∫(integrand_fc)*dΩ_cf
vecdata = collect_cell_vector(U_c,contr_cf_c)
vec_c = assemble_vector(assem_c2f,vecdata)


#end