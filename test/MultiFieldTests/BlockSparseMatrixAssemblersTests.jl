using Gridap
using Gridap.FESpaces, Gridap.Geometry, Gridap.CellData, Gridap.ReferenceFEs, Gridap.Fields

using Gridap.MultiField
using BlockArrays, SparseArrays, LinearAlgebra

sol(x) = sum(x)

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(10,10))
Ω = Triangulation(model)

reffe = LagrangianRefFE(Float64,QUAD,1)
V = FESpace(Ω, reffe; dirichlet_tags="boundary")
U = TrialFESpace(sol,V)

dΩ = Measure(Ω, 2)
biform((u1,u2),(v1,v2)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2)*dΩ
liform((v1,v2)) = ∫(v1 - v2)*dΩ

############################################################################################
# Normal assembly 

Y = MultiFieldFESpace([V,V])
X = MultiFieldFESpace([U,U])

u = get_trial_fe_basis(X)
v = get_fe_basis(Y)

data = collect_cell_matrix_and_vector(X,Y,biform(u,v),liform(v))
matdata = collect_cell_matrix(X,Y,biform(u,v))
vecdata = collect_cell_vector(Y,liform(v))  

mat = assemble_matrix(assem.glob_assembler,matdata)
vec = assemble_vector(assem.glob_assembler,vecdata)

y = similar(vec)
mul!(y,mat,vec)

############################################################################################
# Block Assembly 

mfs = BlockMultiFieldStyle()
Yb = MultiFieldFESpace([V,V];style=mfs)
Xb = MultiFieldFESpace([U,U];style=mfs)

ub = get_trial_fe_basis(Xb)
vb = get_fe_basis(Yb)

bdata = collect_cell_matrix_and_vector(Xb,Yb,biform(ub,vb),liform(vb))
bmatdata = collect_cell_matrix(Xb,Yb,biform(ub,vb))
bvecdata = collect_cell_vector(Yb,liform(vb)) 

assem = BlockSparseMatrixAssembler(Xb,Yb)
mat_blocks = assemble_matrix(assem,bmatdata)
vec_blocks = assemble_vector(assem,bvecdata)

y_blocks = similar(vec_blocks)
mul!(y_blocks,mat_blocks,vec_blocks)

y_blocks ≈ y
