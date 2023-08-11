module BlockMatrixAssemblersTests
using Test, BlockArrays, SparseArrays, LinearAlgebra

using Gridap
using Gridap.FESpaces, Gridap.ReferenceFEs, Gridap.MultiField

function main(n_spaces,mfs,weakform,Ω,dΩ,U,V)
  biform, liform = weakform

  ############################################################################################
  # Normal assembly 

  Y = MultiFieldFESpace(fill(V,n_spaces))
  X = MultiFieldFESpace(fill(U,n_spaces))

  u = get_trial_fe_basis(X)
  v = get_fe_basis(Y)

  data = collect_cell_matrix_and_vector(X,Y,biform(u,v),liform(v))
  matdata = collect_cell_matrix(X,Y,biform(u,v))
  vecdata = collect_cell_vector(Y,liform(v))  

  assem = SparseMatrixAssembler(X,Y)
  A1 = assemble_matrix(assem,matdata)
  b1 = assemble_vector(assem,vecdata)
  A2,b2 = assemble_matrix_and_vector(assem,data)

  ############################################################################################
  # Block MultiFieldStyle

  Yb = MultiFieldFESpace(fill(V,n_spaces);style=mfs)
  Xb = MultiFieldFESpace(fill(U,n_spaces);style=mfs)
  test_fe_space(Yb)
  test_fe_space(Xb)

  ub = get_trial_fe_basis(Xb)
  vb = get_fe_basis(Yb)

  bdata = collect_cell_matrix_and_vector(Xb,Yb,biform(ub,vb),liform(vb))
  bmatdata = collect_cell_matrix(Xb,Yb,biform(ub,vb))
  bvecdata = collect_cell_vector(Yb,liform(vb))
  test_fe_space(Xb,bdata[1][1][1],bmatdata[1][1],bvecdata[1][1],Ω)
  test_fe_space(Yb,bdata[1][1][1],bmatdata[1][1],bvecdata[1][1],Ω)

  ############################################################################################
  # Block Assembly 

  assem_blocks = SparseMatrixAssembler(Xb,Yb)
  test_sparse_matrix_assembler(assem_blocks,bmatdata,bvecdata,bdata)

  A1_blocks = assemble_matrix(assem_blocks,bmatdata)
  b1_blocks = assemble_vector(assem_blocks,bvecdata)
  @test A1 ≈ A1_blocks
  @test b1 ≈ b1_blocks

  y1_blocks = similar(b1_blocks)
  mul!(y1_blocks,A1_blocks,b1_blocks)
  y1 = similar(b1)
  mul!(y1,A1,b1)
  @test y1_blocks ≈ y1

  A2_blocks, b2_blocks = assemble_matrix_and_vector(assem_blocks,bdata)
  @test A2_blocks ≈ A2
  @test b2_blocks ≈ b2

  A3_blocks = allocate_matrix(assem_blocks,bmatdata)
  b3_blocks = allocate_vector(assem_blocks,bvecdata)
  assemble_matrix!(A3_blocks,assem_blocks,bmatdata)
  assemble_vector!(b3_blocks,assem_blocks,bvecdata)
  @test A3_blocks ≈ A1
  @test b3_blocks ≈ b1_blocks

  A4_blocks, b4_blocks = allocate_matrix_and_vector(assem_blocks,bdata)
  assemble_matrix_and_vector!(A4_blocks,b4_blocks,assem_blocks,bdata)
  @test A4_blocks ≈ A2_blocks
  @test b4_blocks ≈ b2_blocks

  ############################################################################################

  op = AffineFEOperator(biform,liform,X,Y)
  block_op = AffineFEOperator(biform,liform,Xb,Yb)

  @test get_matrix(op) ≈ get_matrix(block_op)
  @test get_vector(op) ≈ get_vector(block_op)
end

############################################################################################

sol(x) = sum(x)

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(5,5))
Ω = Triangulation(model)

reffe = LagrangianRefFE(Float64,QUAD,1)
V = FESpace(Ω, reffe; dirichlet_tags="boundary")
U = TrialFESpace(sol,V)

dΩ = Measure(Ω, 2)
biform2((u1,u2),(v1,v2)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2 - u1⋅v2)*dΩ
liform2((v1,v2)) = ∫(v1 - v2)*dΩ
biform3((u1,u2,u3),(v1,v2,v3)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2 - u1⋅v2 - u3⋅v2 - u2⋅v3)*dΩ
liform3((v1,v2,v3)) = ∫(v1 - v2 + 2.0*v3)*dΩ

for (n_spaces,weakform) in zip([2,3],[(biform2,liform2),(biform3,liform3)])
  for mfs in [BlockMultiFieldStyle(),BlockMultiFieldStyle(2,(1,n_spaces-1))]
    main(n_spaces,mfs,weakform,Ω,dΩ,U,V)
  end
end

end # module