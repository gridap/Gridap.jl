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
using Gridap.FESpaces
using Gridap.CellData

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

reffe = ReferenceFE(:Lagrangian,Float64,1)
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

bcellmat = integrate(v*u,bquad)
bcellvec = integrate(v*3,bquad)
bcellmatvec = pair_arrays(bcellmat,bcellvec)
bcellids = btrian.glue.face_to_cell

b0cellmat = integrate(v*u,b0quad)
b0cellvec = integrate(v*3,b0quad)
b0cellmatvec = pair_arrays(b0cellmat,b0cellvec)
b0cellids = get_cell_id(b0trian)
@test length(b0cellids) == 0

term_to_cellmat = [cellmat, bcellmat, b0cellmat]
term_to_cellvec = [cellvec, bcellvec, b0cellvec]
term_to_cellids = [cellids, bcellids, b0cellids]
term_to_cellmatvec = [ cellmatvec, bcellmatvec, b0cellmatvec ]

mtypes = [
  SparseMatrixCSC,
  SparseMatrixCSR,
  SparseMatrixCSR{0},
  SparseMatrixCSR{1},
  SparseMatrixCSR{0,Float64, Int},
  SparseMatrixCSR{1,Float64, Int},
  SymSparseMatrixCSR,
  SymSparseMatrixCSR{0},
  SymSparseMatrixCSR{1},
  SymSparseMatrixCSR{0, Float64, Int},
  SymSparseMatrixCSR{1, Float64, Int}]

for T in mtypes

  matvecdata = ( term_to_cellmatvec , term_to_cellids, term_to_cellids)
  matdata = (term_to_cellmat,term_to_cellids,term_to_cellids)
  vecdata = (term_to_cellvec,term_to_cellids)
  data = (matvecdata,matdata,vecdata)

  assem = SparseMatrixAssembler(T,Vector{Float64},U,V)
  test_sparse_matrix_assembler(assem,matdata,vecdata,data)
  
  matdata = ([cellmat],[cellids],[cellids])
  vecdata = ([cellvec],[cellids])

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

  data = (([cellmatvec],[cellids],[cellids]),([],[],[]),([],[]))
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
scellids = get_cell_id(strian)
zh = zero(V)
scellvals = get_cell_dof_values(zh,scellids)
scellmatvec = attach_dirichlet(scellmatvec,scellvals)

matvecdata = ([scellmatvec],[scellids],[scellids])
matdata = ([scellmat],[scellids],[scellids])
vecdata = ([scellvec],[scellids])
data = (matvecdata,matdata,vecdata)

assem = SparseMatrixAssembler(U,V)
test_sparse_matrix_assembler(assem,matdata,vecdata,data)

A = assemble_matrix(assem,matdata)
@test A == zeros(num_free_dofs(V),num_free_dofs(U))

# Allowing a non-blocked matrix in a context where you would expect a blocked one.
scellmat = map(i->zeros(size(i)),scellmat)
matdata = ([scellmat],[scellids],[scellids])
A = assemble_matrix(assem,matdata)
@test A == zeros(num_free_dofs(V),num_free_dofs(U))

end # module
