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

assem = SparseMatrixAssembler(X,Y)
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)

y = similar(b)
mul!(y,A,b)

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

assem_blocks = SparseMatrixAssembler(Xb,Yb)
A_blocks = assemble_matrix(assem_blocks,bmatdata)
b_blocks = assemble_vector(assem_blocks,bvecdata)

y_blocks = similar(b_blocks)
mul!(y_blocks,A_blocks,b_blocks)

y_blocks ≈ y
