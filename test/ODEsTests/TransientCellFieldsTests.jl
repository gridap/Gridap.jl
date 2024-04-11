module TransientCellFieldsTests

using Test
using LinearAlgebra
using SparseArrays
using BlockArrays

using Gridap
using Gridap.CellData
using Gridap.FESpaces
using Gridap.MultiField
using Gridap.ODEs
using Gridap.ODEs: TransientMultiFieldCellField

f(t) = x -> sum(x)

domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, f)

Ω = Triangulation(model)
degree = 2 * order
dΩ = Measure(Ω, degree)

###############
# SingleField #
###############
m(t, ∂ₜu, v) = ∫(∂ₜu ⋅ v) * dΩ
b(t, u, v) = ∫(∇(u) ⋅ ∇(v) + u ⋅ v) * dΩ
l(t, v) = ∫(v) * dΩ

res(t, u, v) = m(t, ∂t(u), v) + b(t, u, v) - l(t, v)
jac(t, u, du, v) = b(t, du, v)
jac_t(t, u, dut, v) = m(t, dut, v)

t0 = 0.0
U0 = U(t0)
u = get_trial_fe_basis(U0)
uₜ = TransientCellField(u, (u,))
v = get_fe_basis(V)

dc = DomainContribution()
dc = dc + jac(t0, uₜ, u, v)
dc = dc + jac_t(t0, uₜ, u, v)
matdata = collect_cell_matrix(U0, V, dc)
vecdata = collect_cell_vector(V, l(t0, v))

assembler = SparseMatrixAssembler(U, V)
mat = assemble_matrix(assembler, matdata)
vec = assemble_vector(assembler, vecdata)

##############
# MultiField #
##############
m2(t, (∂ₜu1, ∂ₜu2), (v1, v2)) = ∫(∂ₜu1 ⋅ v1) * dΩ
b2(t, (u1, u2), (v1, v2)) = ∫(∇(u1) ⋅ ∇(v1) + u2 ⋅ v2 - u1 ⋅ v2) * dΩ
l2(t, (v1, v2)) = ∫(v1 - v2) * dΩ

m3(t, (∂ₜu1, ∂ₜu2, ∂ₜu3), (v1, v2, v3)) = ∫(∂ₜu1 ⋅ v1) * dΩ
b3(t, (u1, u2, u3), (v1, v2, v3)) = ∫(∇(u1) ⋅ ∇(v1) + u2 ⋅ v2 - u1 ⋅ v2 - u3 ⋅ v2 - u2 ⋅ v3) * dΩ
l3(t, (v1, v2, v3)) = ∫(v1 - v2 + v3) * dΩ

function test_multifield(n, mfs, m, b, l, U, V)
  res(t, u, v) = m(t, ∂t(u), v) + b(t, u, v) - l(t, v)
  jac(t, u, du, v) = b(t, du, v)
  jac_t(t, u, dut, v) = m(t, dut, v)

  # Normal assembly
  Y = MultiFieldFESpace(fill(V, n))
  X = TransientMultiFieldFESpace(fill(U, n))

  t0 = 0.0
  X0 = X(t0)
  test_fe_space(Y)
  test_fe_space(X0)

  u = get_trial_fe_basis(X0)
  uₜ = TransientMultiFieldCellField(u, (u,))
  v = get_fe_basis(Y)

  dc = DomainContribution()
  dc = dc + jac(t0, uₜ, u, v)
  dc = dc + jac_t(t0, uₜ, u, v)
  matdata = collect_cell_matrix(X0, Y, dc)
  vecdata = collect_cell_vector(Y, l(t0, v))

  assembler = SparseMatrixAssembler(X, Y)
  mat = assemble_matrix(assembler, matdata)
  vec = assemble_vector(assembler, vecdata)

  # Block MultiFieldStyle
  Y_blocks = MultiFieldFESpace(fill(V, n); style=mfs)
  X_blocks = TransientMultiFieldFESpace(fill(U, n); style=mfs)
  X0_blocks = X_blocks(t0)
  test_fe_space(Y_blocks)
  test_fe_space(X0_blocks)

  u_blocks = get_trial_fe_basis(X0_blocks)
  uₜ_blocks = TransientMultiFieldCellField(u_blocks, (u_blocks,))
  v_blocks = get_fe_basis(Y_blocks)

  dc = DomainContribution()
  dc = dc + jac(t0, uₜ_blocks, u_blocks, v_blocks)
  dc = dc + jac_t(t0, uₜ_blocks, u_blocks, v_blocks)
  matdata_blocks = collect_cell_matrix(X0_blocks, Y_blocks, dc)
  vecdata_blocks = collect_cell_vector(Y_blocks, l(t0, v_blocks))

  # Block Assembly
  assembler_blocks = SparseMatrixAssembler(X_blocks, Y_blocks)

  mat_blocks = assemble_matrix(assembler_blocks, matdata_blocks)
  vec_blocks = assemble_vector(assembler_blocks, vecdata_blocks)
  @test mat_blocks ≈ mat
  @test vec_blocks ≈ vec

  matvec = similar(vec)
  mul!(matvec, mat, vec)
  matvec_blocks = similar(vec_blocks)
  mul!(matvec_blocks, mat_blocks, vec_blocks)
  @test matvec_blocks ≈ matvec

  mat_blocks = allocate_matrix(assembler_blocks, matdata_blocks)
  vec_blocks = allocate_vector(assembler_blocks, vecdata_blocks)
  assemble_matrix!(mat_blocks, assembler_blocks, matdata_blocks)
  assemble_vector!(vec_blocks, assembler_blocks, vecdata_blocks)
  @test mat_blocks ≈ mat
  @test vec_blocks ≈ vec
end

for (n, m, b, l) in ((2, m2, b2, l2), (3, m3, b3, l3),)
  for mfs in (BlockMultiFieldStyle(), BlockMultiFieldStyle(2, (1, n - 1)))
    test_multifield(n, mfs, m, b, l, U, V)
  end
end

end # module TransientCellFieldsTests
