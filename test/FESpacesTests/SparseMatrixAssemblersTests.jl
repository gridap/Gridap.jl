module SparseMatrixAssemblers

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.Algebra
using SparseArrays
using SparseMatricesCSR
using Gridap.FESpaces
using Gridap.CellData

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(model,reffe,dirichlet_tags=[1,2,3,4,6,5])
U = V

b(x) = x[2]

v = get_cell_shapefuns(V)
u = get_cell_shapefuns_trial(U)

degree = 2
trian = get_triangulation(model)
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,degree)

b0trian = Triangulation(btrian,Int[])
b0quad = CellQuadrature(b0trian,degree)

cellmat = integrate(∇(v)⊙∇(u),quad)
cellvec = integrate(v⊙b,quad)
cellmatvec = pair_arrays(cellmat,cellvec)
cellids = collect(1:num_cells(trian))
rows = get_cell_dof_ids(V,cellids)
cols = get_cell_dof_ids(U,cellids)
cellmat_c = attach_constraints_cols(U,cellmat,cellids)
cellmat_rc = attach_constraints_rows(V,cellmat_c,cellids)
cellvec_r = attach_constraints_rows(V,cellvec,cellids)
cellmatvec_c = attach_constraints_cols(U,cellmatvec,cellids)
cellmatvec_rc = attach_constraints_rows(V,cellmatvec_c,cellids)

bcellmat = integrate(v*u,bquad)
bcellvec = integrate(v*3,bquad)
bcellmatvec = pair_arrays(bcellmat,bcellvec)
bcellids = btrian.glue.face_to_cell
brows = get_cell_dof_ids(V,bcellids)
bcols = get_cell_dof_ids(U,bcellids)
bcellmat_c = attach_constraints_cols(U,bcellmat,bcellids)
bcellmat_rc = attach_constraints_rows(V,bcellmat_c,bcellids)
bcellvec_r = attach_constraints_rows(V,bcellvec,bcellids)
bcellmatvec_c = attach_constraints_cols(U,bcellmatvec,bcellids)
bcellmatvec_rc = attach_constraints_rows(V,bcellmatvec_c,bcellids)

b0cellmat = integrate(v*u,b0quad)
b0cellvec = integrate(v*3,b0quad)
b0cellmatvec = pair_arrays(b0cellmat,b0cellvec)
b0cellids = get_cell_to_bgcell(b0trian)
b0rows = get_cell_dof_ids(V,b0cellids)
b0cols = get_cell_dof_ids(U,b0cellids)
@test length(b0cellids) == 0
b0cellmat_c = attach_constraints_cols(U,b0cellmat,b0cellids)
b0cellmat_rc = attach_constraints_rows(V,b0cellmat_c,b0cellids)
b0cellvec_r = attach_constraints_rows(V,b0cellvec,b0cellids)
b0cellmatvec_c = attach_constraints_cols(U,b0cellmatvec,b0cellids)
b0cellmatvec_rc = attach_constraints_rows(V,b0cellmatvec_c,b0cellids)

term_to_cellmat = [cellmat_rc, bcellmat_rc, b0cellmat_rc]
term_to_cellvec = [cellvec, bcellvec, b0cellvec]
term_to_rows = [rows, brows, b0rows]
term_to_cols = [cols, bcols, b0cols]
term_to_cellmatvec = [ cellmatvec, bcellmatvec, b0cellmatvec ]

struct AssemblyStrategyMock <: AssemblyStrategy end
FESpaces.row_map(a::AssemblyStrategyMock,row) = row
FESpaces.col_map(a::AssemblyStrategyMock,col) = col
FESpaces.row_mask(a::AssemblyStrategyMock,row) = true
FESpaces.col_mask(a::AssemblyStrategyMock,col) = true

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

  assem2 = SparseMatrixAssembler(T,Vector{Float64},U,V,AssemblyStrategyMock())
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
@test isa(scellmat[1],BlockArrayCoo)
@test is_zero_block(scellmat[1],1,1)
@test is_zero_block(scellmat[1],2,1)
@test is_nonzero_block(scellmat[1],1,2)
@test is_nonzero_block(scellmat[1],2,2)

scellvec = integrate(mean(v*3),squad)
@test isa(scellvec[1],BlockArrayCoo)
scellmatvec = pair_arrays(scellmat,scellvec)
scellids = get_cell_to_bgcell(strian)
srows = get_cell_dof_ids(V,scellids)
scols = get_cell_dof_ids(U,scellids)
zh = zero(V)
scellvals = get_cell_dof_values(zh,scellids)
scellmatvec = attach_dirichlet(scellmatvec,scellvals)

matvecdata = ([scellmatvec],[srows],[scols])
matdata = ([scellmat],[srows],[scols])
vecdata = ([scellvec],[srows])
data = (matvecdata,matdata,vecdata)

assem = SparseMatrixAssembler(U,V)
test_sparse_matrix_assembler(assem,matdata,vecdata,data)

A = assemble_matrix(assem,matdata)
@test A == zeros(num_free_dofs(V),num_free_dofs(U))

# Allowing a non-blocked matrix in a context where you would expect a blocked one.
scellmat = map(i->zeros(size(i)),scellmat)
matdata = ([scellmat],[srows],[scols])
A = assemble_matrix(assem,matdata)
@test A == zeros(num_free_dofs(V),num_free_dofs(U))

end # module
