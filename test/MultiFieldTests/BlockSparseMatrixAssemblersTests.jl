using Gridap
using Gridap.FESpaces, Gridap.Geometry, Gridap.CellData, Gridap.ReferenceFEs, Gridap.Fields

using Gridap.MultiField
using BlockArrays, SparseArrays, LinearAlgebra

############################################################################################

############################################################################################

sol(x) = sum(x)

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(10,10))
Ω = Triangulation(model)

reffe = LagrangianRefFE(Float64,QUAD,1)
V = FESpace(Ω, reffe; dirichlet_tags="boundary")
U = TrialFESpace(sol,V)

Y = MultiFieldFESpace([V,V])
X = MultiFieldFESpace([U,U])

dΩ = Measure(Ω, 2)
biform((u1,u2),(v1,v2)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2)*dΩ
liform((v1,v2)) = ∫(v1 - v2)*dΩ

op = AffineFEOperator(biform,liform,X,Y)

u = get_trial_fe_basis(X)
v = get_fe_basis(Y)

data = collect_cell_matrix_and_vector(X,Y,biform(u,v),liform(v))
matdata = collect_cell_matrix(X,Y,biform(u,v))
vecdata = collect_cell_vector(Y,liform(v))  

############################################################################################
# Block Assembly 
assem = BlockSparseMatrixAssembler(X,Y)
mat_blocks = assemble_matrix(assem,matdata)
vec_blocks = assemble_vector(assem,vecdata)

y_blocks = similar(vec_blocks)
mul!(y_blocks,mat_blocks,vec_blocks)

# Normal assembly 
mat = assemble_matrix(assem.glob_assembler,matdata)
vec = assemble_vector(assem.glob_assembler,vecdata)

y = similar(vec)
mul!(y,mat,vec)

y_blocks ≈ y
