module TransientBlockMultiFieldStyleTests
using Test, BlockArrays, SparseArrays, LinearAlgebra

using Gridap
using Gridap.FESpaces, Gridap.ReferenceFEs, Gridap.MultiField
using Gridap.ODEs.TransientFETools

function main(n_spaces,mfs,weakform,Ω,dΩ,U,V)
  mass, biform, liform = weakform
  res(t,x,y) = mass(t,∂t(x),y) + biform(t,x,y) - liform(t,y)
  jac(t,x,dx,y) = biform(t,dx,y)
  jac_t(t,xt,dxt,y) = mass(t,dxt,y)

  ############################################################################################
  # Normal assembly

  Y = MultiFieldFESpace(fill(V,n_spaces))
  X = TransientMultiFieldFESpace(fill(U,n_spaces))

  u = get_trial_fe_basis(X(0.0))
  v = get_fe_basis(Y)
  uₜ = TransientCellField(u,(u,))

  matdata_jac = collect_cell_matrix(X(0),Y,jac(0,uₜ,u,v))
  matdata_jac_t = collect_cell_matrix(X(0),Y,jac_t(0,uₜ,u,v))
  matdata_jacs = (matdata_jac,matdata_jac_t)
  matdata = TransientFETools._vcat_matdata(matdata_jacs)
  vecdata = collect_cell_vector(Y,liform(0,v))

  assem = SparseMatrixAssembler(X(0),Y)
  A1 = assemble_matrix(assem,matdata)
  b1 = assemble_vector(assem,vecdata)

  ############################################################################################
  # Block MultiFieldStyle

  Yb = MultiFieldFESpace(fill(V,n_spaces);style=mfs)
  Xb = TransientMultiFieldFESpace(fill(U,n_spaces);style=mfs)
  test_fe_space(Yb)
  test_fe_space(Xb(0))

  ub = get_trial_fe_basis(Xb(0))
  vb = get_fe_basis(Yb)
  ubₜ = TransientCellField(ub,(ub,))

  bmatdata_jac = collect_cell_matrix(Xb(0),Yb,jac(0,ubₜ,ub,vb))
  bmatdata_jac_t = collect_cell_matrix(Xb(0),Yb,jac_t(0,ubₜ,ub,vb))
  bmatdata_jacs = (bmatdata_jac,bmatdata_jac_t)
  bmatdata = TransientFETools._vcat_matdata(bmatdata_jacs)
  bvecdata = collect_cell_vector(Yb,liform(0,vb))

  ############################################################################################
  # Block Assembly

  assem_blocks = SparseMatrixAssembler(Xb,Yb)

  A1_blocks = assemble_matrix(assem_blocks,bmatdata)
  b1_blocks = assemble_vector(assem_blocks,bvecdata)
  @test A1 ≈ A1_blocks
  @test b1 ≈ b1_blocks

  y1_blocks = similar(b1_blocks)
  mul!(y1_blocks,A1_blocks,b1_blocks)
  y1 = similar(b1)
  mul!(y1,A1,b1)
  @test y1_blocks ≈ y1

  A3_blocks = allocate_matrix(assem_blocks,bmatdata)
  b3_blocks = allocate_vector(assem_blocks,bvecdata)
  assemble_matrix!(A3_blocks,assem_blocks,bmatdata)
  assemble_vector!(b3_blocks,assem_blocks,bvecdata)
  @test A3_blocks ≈ A1
  @test b3_blocks ≈ b1_blocks

end

############################################################################################

sol(x,t) = sum(x)
sol(t::Real) = x->sol(x,t)

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(5,5))
Ω = Triangulation(model)

reffe = LagrangianRefFE(Float64,QUAD,1)
V = FESpace(Ω, reffe; dirichlet_tags="boundary")
U = TransientTrialFESpace(V,sol)

dΩ = Measure(Ω, 2)
mass2(t,(u1t,u2t),(v1,v2)) = ∫(u1t⋅v1)*dΩ
biform2(t,(u1,u2),(v1,v2)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2 - u1⋅v2)*dΩ
liform2(t,(v1,v2)) = ∫(v1 - v2)*dΩ
mass3(t,(u1t,u2t,u3t),(v1,v2,v3)) = ∫(u1t⋅v1)*dΩ
biform3(t,(u1,u2,u3),(v1,v2,v3)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2 - u1⋅v2 - u3⋅v2 - u2⋅v3)*dΩ
liform3(t,(v1,v2,v3)) = ∫(v1 - v2 + 2.0*v3)*dΩ

for (n_spaces,weakform) in zip([2,3],[(mass2,biform2,liform2),(mass3,biform3,liform3)])
  for mfs in [BlockMultiFieldStyle(),BlockMultiFieldStyle(2,(1,n_spaces-1))]
    main(n_spaces,mfs,weakform,Ω,dΩ,U,V)
  end
end

end # module
