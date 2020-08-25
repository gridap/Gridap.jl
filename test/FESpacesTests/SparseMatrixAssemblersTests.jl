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

order = 1
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

dirichlet_tags = [1,2,3,4,6,5]
V = GradConformingFESpace(reffes,model,dirichlet_tags)

U = TrialFESpace(V)

b(x) = x[2]

v = get_cell_basis(V)
u = get_cell_basis(U)

trian = get_triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,degree)

bu = restrict(u,btrian)
bv = restrict(v,btrian)

cellmat = integrate(∇(v)⊙∇(u),trian,quad)
cellvec = integrate(v⊙b,trian,quad)
cellmatvec = pair_arrays(cellmat,cellvec)
cellids = collect(1:num_cells(trian))

bcellmat = integrate(bv*bu,btrian,bquad)
bcellvec = integrate(bv*3,btrian,bquad)
bcellmatvec = pair_arrays(bcellmat,bcellvec)
bcellids = get_face_to_cell(btrian)

term_to_cellmat = [cellmat, bcellmat]
term_to_cellvec = [cellvec, bcellvec]
term_to_cellids = [cellids, bcellids]
term_to_cellmatvec = [ cellmatvec, bcellmatvec ]


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

su = restrict(u,strian)
sv = restrict(v,strian)

scellmat = integrate(jump(sv)*su.⁻,strian,squad)
scellvec = integrate(mean(sv*3),strian,squad)
scellmatvec = pair_arrays(scellmat,scellvec)
scellids = get_cell_id(strian)
zh = zero(V)
scellvals = get_cell_values(zh,scellids)
scellmatvec = attach_dirichlet(scellmatvec,scellvals)

matvecdata = ([scellmatvec],[scellids],[scellids])
matdata = ([scellmat],[scellids],[scellids])
vecdata = ([scellvec],[scellids])
data = (matvecdata,matdata,vecdata)

assem = SparseMatrixAssembler(U,V)
test_sparse_matrix_assembler(assem,matdata,vecdata,data)

A = assemble_matrix(assem,matdata)
@test A == zeros(num_free_dofs(V),num_free_dofs(U))

end # module
