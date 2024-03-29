module SparseMatrixAssemblers

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.Algebra
using SparseArrays
using SparseMatricesCSR
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Algebra

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(model,reffe,dirichlet_tags=[1,2,3,4,6,5])
U = V

b(x) = x[2]

v = get_fe_basis(V)
u = get_trial_fe_basis(U)

degree = 2
trian = get_triangulation(model)
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,degree)

b0trian = Triangulation(model,Int[])
b0quad = CellQuadrature(b0trian,degree)

cellmat = integrate(∇(v)⊙∇(u),quad)
cellvec = integrate(v⊙b,quad)
cellmatvec = pair_arrays(cellmat,cellvec)
rows = get_cell_dof_ids(V,trian)
cols = get_cell_dof_ids(U,trian)
cellmat_c = attach_constraints_cols(U,cellmat,trian)
cellmat_rc = attach_constraints_rows(V,cellmat_c,trian)
cellvec_r = attach_constraints_rows(V,cellvec,trian)
cellmatvec_c = attach_constraints_cols(U,cellmatvec,trian)
cellmatvec_rc = attach_constraints_rows(V,cellmatvec_c,trian)

bcellmat = integrate(v*u,bquad)
bcellvec = integrate(v*3,bquad)
bcellmatvec = pair_arrays(bcellmat,bcellvec)
brows = get_cell_dof_ids(V,btrian)
bcols = get_cell_dof_ids(U,btrian)
bcellmat_c = attach_constraints_cols(U,bcellmat,btrian)
bcellmat_rc = attach_constraints_rows(V,bcellmat_c,btrian)
bcellvec_r = attach_constraints_rows(V,bcellvec,btrian)
bcellmatvec_c = attach_constraints_cols(U,bcellmatvec,btrian)
bcellmatvec_rc = attach_constraints_rows(V,bcellmatvec_c,btrian)

b0cellmat = integrate(v*u,b0quad)
b0cellvec = integrate(v*3,b0quad)
b0cellmatvec = pair_arrays(b0cellmat,b0cellvec)
b0rows = get_cell_dof_ids(V,b0trian)
b0cols = get_cell_dof_ids(U,b0trian)
@test length(b0rows) == 0
b0cellmat_c = attach_constraints_cols(U,b0cellmat,b0trian)
b0cellmat_rc = attach_constraints_rows(V,b0cellmat_c,b0trian)
b0cellvec_r = attach_constraints_rows(V,b0cellvec,b0trian)
b0cellmatvec_c = attach_constraints_cols(U,b0cellmatvec,b0trian)
b0cellmatvec_rc = attach_constraints_rows(V,b0cellmatvec_c,b0trian)

term_to_cellmat = [cellmat_rc, bcellmat_rc, b0cellmat_rc]
term_to_cellvec = [cellvec, bcellvec, b0cellvec]
term_to_rows = [rows, brows, b0rows]
term_to_cols = [cols, bcols, b0cols]
term_to_cellmatvec = [ cellmatvec, bcellmatvec, b0cellmatvec ]

mtypes = [
  SparseMatrixCSC{Float64,Int},
  SparseMatrixCSR{0,Float64,Int},
  SparseMatrixCSR{1,Float64,Int},
  SymSparseMatrixCSR{0,Float64,Int},
  SymSparseMatrixCSR{1,Float64,Int}]

for T in mtypes

  matvecdata = ( term_to_cellmatvec , term_to_rows, term_to_cols)
  matdata = (term_to_cellmat,term_to_rows,term_to_cols)
  vecdata = (term_to_cellvec,term_to_rows)
  data = (matvecdata,matdata,vecdata)

  assem = SparseMatrixAssembler(T,Vector{Float64},U,V)
  test_sparse_matrix_assembler(assem,matdata,vecdata,data)

  strategy = GenericAssemblyStrategy(row->row,col->col,row->true,col->true)

  assem2 = SparseMatrixAssembler(T,Vector{Float64},U,V,strategy)
  test_sparse_matrix_assembler(assem2,matdata,vecdata,data)

  matdata = ([cellmat],[rows],[cols])
  vecdata = ([cellvec],[rows])

  mat = assemble_matrix(assem,matdata)
  vec = assemble_vector(assem,vecdata)

  x = mat \ vec

  assemble_matrix!(mat,assem,matdata)
  assemble_vector!(vec,assem,vecdata)

  x2 = mat \ vec
  @test x ≈ x2

  @test vec ≈ [0.0625, 0.125, 0.0625]
  @test mat[1, 1]  ≈  1.333333333333333
  @test mat[2, 1]  ≈ -0.33333333333333
  @test mat[1, 2]  ≈ -0.33333333333333
  @test mat[2, 2]  ≈ 2.666666666666666
  @test mat[3, 2]  ≈ -0.33333333333333
  @test mat[2, 3]  ≈ -0.33333333333333
  @test mat[3, 3]  ≈ 1.333333333333333

  data = (([cellmatvec],[rows],[cols]),([],[],[]),([],[]))
  mat, vec = allocate_matrix_and_vector(assem,data)
  assemble_matrix_and_vector!(mat,vec,assem,data)
  assemble_matrix_and_vector!(mat,vec,assem,data)

  @test vec ≈ [0.0625, 0.125, 0.0625]
  @test mat[1, 1]  ≈  1.333333333333333
  @test mat[2, 1]  ≈ -0.33333333333333

  x3 = mat \ vec
  @test x ≈ x3

  mat, vec = assemble_matrix_and_vector(assem,data)

  x4 = mat \ vec
  @test x ≈ x4

  @test vec ≈ [0.0625, 0.125, 0.0625]
  @test mat[1, 1]  ≈  1.333333333333333
  @test mat[2, 1]  ≈ -0.33333333333333

  mat, vec = assemble_matrix_and_vector(assem2,data)

  x4 = mat \ vec
  @test x ≈ x4

  @test vec ≈ [0.0625, 0.125, 0.0625]
  @test mat[1, 1]  ≈  1.333333333333333
  @test mat[2, 1]  ≈ -0.33333333333333

end

strian = SkeletonTriangulation(model)
squad = CellQuadrature(strian,degree)

scellmat = integrate(jump(v)*u.⁻,squad)
@test isa(scellmat[1],ArrayBlock)
@test scellmat[1][1,1] === nothing
@test scellmat[1][2,1] === nothing
@test scellmat[1][1,2] !== nothing
@test scellmat[1][2,2] !== nothing

scellvec = integrate(mean(v*3),squad)
@test isa(scellvec[1],ArrayBlock)
scellmatvec = pair_arrays(scellmat,scellvec)
srows = get_cell_dof_ids(V,strian)
scols = get_cell_dof_ids(U,strian)
zh = zero(V)
scellvals = get_cell_dof_values(zh,strian)
scellmatvec = attach_dirichlet(scellmatvec,scellvals)

matvecdata = ([scellmatvec],[srows],[scols])
matdata = ([scellmat],[srows],[scols])
vecdata = ([scellvec],[srows])
data = (matvecdata,matdata,vecdata)

assem = SparseMatrixAssembler(U,V)
test_sparse_matrix_assembler(assem,matdata,vecdata,data)

A = assemble_matrix(assem,matdata)
@test A == zeros(num_free_dofs(V),num_free_dofs(U))

assem = SparseMatrixAssembler(
  SparseMatrixBuilder(SparseMatrixCSC{Float64,Int},MinCPU()),U,V)
test_sparse_matrix_assembler(assem,matdata,vecdata,data)

assem = SparseMatrixAssembler(
  SparseMatrixBuilder(SparseMatrixCSC{Float64,Int},MinMemory()),U,V)
test_sparse_matrix_assembler(assem,matdata,vecdata,data)

end # module
