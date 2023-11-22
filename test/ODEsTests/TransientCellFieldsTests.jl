module TransientCellFieldsTests

using Test
using LinearAlgebra
using SparseArrays
using BlockArrays

using Gridap
using Gridap.ODEs

u(x, t) = sum(x)
u(t::Real) = x -> u(x, t)

domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, u)

Ω = Triangulation(model)
degree = 2
dΩ = Measure(Ω, degree)

m1(t, u1t, v1) = ∫(u1t ⋅ v1) * dΩ
b1(t, u1, v1) = ∫(∇(u1) ⋅ ∇(v1) + u1 ⋅ v1) * dΩ
l1(t, v1) = ∫(v1) * dΩ

m2(t, (u1t, u2t), (v1, v2)) = ∫(u1t ⋅ v1) * dΩ
b2(t, (u1, u2), (v1, v2)) = ∫(∇(u1) ⋅ ∇(v1) + u2 ⋅ v2 - u1 ⋅ v2) * dΩ
l2(t, (v1, v2)) = ∫(v1 - v2) * dΩ

m3(t, (u1t, u2t, u3t), (v1, v2, v3)) = ∫(u1t ⋅ v1) * dΩ
b3(t, (u1, u2, u3), (v1, v2, v3)) = ∫(∇(u1) ⋅ ∇(v1) + u2 ⋅ v2 - u1 ⋅ v2 - u3 ⋅ v2 - u2 ⋅ v3) * dΩ
l3(t, (v1, v2, v3)) = ∫(v1 - v2 + v3) * dΩ

# (1, m1, b1, l1)

function test_multifield(n, mfs, m, b, l, Ω, dΩ, U, V)
  mass, biform, liform = weakform
  res(t, x, y) = mass(t, ∂t(x), y) + biform(t, x, y) - liform(t, y)
  jac(t, x, dx, y) = biform(t, dx, y)
  jac_t(t, xt, dxt, y) = mass(t, dxt, y)

  # Normal assembly
  Y = MultiFieldFESpace(fill(V, n_spaces))
  X = TransientMultiFieldFESpace(fill(U, n_spaces))

  u = get_trial_fe_basis(X(0.0))
  v = get_fe_basis(Y)
  uₜ = TransientCellField(u, (u,))

  matdata_jac = collect_cell_matrix(X(0), Y, jac(0, uₜ, u, v))
  matdata_jac_t = collect_cell_matrix(X(0), Y, jac_t(0, uₜ, u, v))
  matdata_jacs = (matdata_jac, matdata_jac_t)
  matdata = TransientFETools._vcat_matdata(matdata_jacs)
  vecdata = collect_cell_vector(Y, liform(0, v))

  assem = SparseMatrixAssembler(X(0), Y)
  A1 = assemble_matrix(assem, matdata)
  b1 = assemble_vector(assem, vecdata)

  # Block MultiFieldStyle
  Yb = MultiFieldFESpace(fill(V, n_spaces); style=mfs)
  Xb = TransientMultiFieldFESpace(fill(U, n_spaces); style=mfs)
  test_fe_space(Yb)
  test_fe_space(Xb(0))

  ub = get_trial_fe_basis(Xb(0))
  vb = get_fe_basis(Yb)
  ubₜ = TransientCellField(ub, (ub,))

  bmatdata_jac = collect_cell_matrix(Xb(0), Yb, jac(0, ubₜ, ub, vb))
  bmatdata_jac_t = collect_cell_matrix(Xb(0), Yb, jac_t(0, ubₜ, ub, vb))
  bmatdata_jacs = (bmatdata_jac, bmatdata_jac_t)
  bmatdata = TransientFETools._vcat_matdata(bmatdata_jacs)
  bvecdata = collect_cell_vector(Yb, liform(0, vb))

  # Block Assembly
  assem_blocks = SparseMatrixAssembler(Xb, Yb)

  A1_blocks = assemble_matrix(assem_blocks, bmatdata)
  b1_blocks = assemble_vector(assem_blocks, bvecdata)
  @test A1 ≈ A1_blocks
  @test b1 ≈ b1_blocks

  y1_blocks = similar(b1_blocks)
  mul!(y1_blocks, A1_blocks, b1_blocks)
  y1 = similar(b1)
  mul!(y1, A1, b1)
  @test y1_blocks ≈ y1

  A3_blocks = allocate_matrix(assem_blocks, bmatdata)
  b3_blocks = allocate_vector(assem_blocks, bvecdata)
  assemble_matrix!(A3_blocks, assem_blocks, bmatdata)
  assemble_vector!(b3_blocks, assem_blocks, bvecdata)
  @test A3_blocks ≈ A1
  @test b3_blocks ≈ b1_blocks
end

for (n, m, b, l) in ((2, m2, b2, l2), (3, m3, b3, l3),)
  for mfs in (BlockMultiFieldStyle(), BlockMultiFieldStyle(2, (1, n - 1)))
    test_multifield(n, mfs, m, b, l, Ω, dΩ, U, V)
  end
end

end # module TransientCellFieldsTests
