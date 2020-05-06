
"""
"""
abstract type AssemblyStrategy end

"""
"""
function row_map(a::AssemblyStrategy,row)
  @abstractmethod
end

"""
"""
function col_map(a::AssemblyStrategy,col)
  @abstractmethod
end

"""
"""
function row_mask(a::AssemblyStrategy,row)
  @abstractmethod
end

"""
"""
function col_mask(a::AssemblyStrategy,col)
  @abstractmethod
end

struct DefaultAssemblyStrategy <: AssemblyStrategy end

row_map(a::DefaultAssemblyStrategy,row) = row

col_map(a::DefaultAssemblyStrategy,col) = col

row_mask(a::DefaultAssemblyStrategy,row) = true

col_mask(a::DefaultAssemblyStrategy,col) = true

"""
"""
abstract type Assembler <: GridapType end

"""
"""
function get_test(a::Assembler)
  @abstractmethod
end

"""
"""
function get_trial(a::Assembler)
  @abstractmethod
end

"""
"""
function get_assembly_strategy(a::Assembler)
  @abstractmethod
end

"""
"""
function allocate_matrix(a::Assembler,cellidsrows,cellidscols)
  @abstractmethod
end

function allocate_matrix(a::Assembler,cellmat,cellidsrows,cellidscols)
  allocate_matrix(a,cellidsrows,cellidscols)
end

"""
"""
function allocate_vector(a::Assembler,cellidsrows)
  @abstractmethod
end

function allocate_vector(a::Assembler,cellvec,cellidsrows)
  allocate_vector(a,cellidsrows)
end

"""
"""
function allocate_matrix_and_vector(a::Assembler,matvecdata,matdata,vecdata)
  matrows, matcols, vecrows = _rearange_cell_ids(matvecdata,matdata,vecdata)
  A = allocate_matrix(a,matrows,matcols)
  b = allocate_vector(a,vecrows)
  (A,b)
end

function allocate_matrix_and_vector(a::Assembler,matvecdata)
  matdata = ([],[],[])
  vecdata = ([],[])
  allocate_matrix_and_vector(a,matvecdata,matdata,vecdata)
end

"""
"""
function assemble_matrix!(A,a::Assembler,cellmat,cellidsrows,cellidscols)
  @abstractmethod
end

"""
"""
function assemble_matrix_add!(
  mat,a::Assembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
  @abstractmethod
end

"""
"""
function assemble_vector!(b,a::Assembler,cellvec,cellids)
  @abstractmethod
end

"""
"""
function assemble_vector_add!(b,a::Assembler,term_to_cellvec,term_to_cellidsrows)
  @abstractmethod
end

"""
"""
function assemble_matrix_and_vector!(A,b,a::Assembler,matvecdata,matdata,vecdata)
  @abstractmethod
end

function assemble_matrix_and_vector_add!(A,b,a::Assembler, matvecdata, matdata, vecdata)
  @abstractmethod
end

function assemble_matrix_and_vector!(A,b,a::Assembler,matvecdata)
  matdata = ([],[],[])
  vecdata = ([],[])
  assemble_matrix_and_vector!(A,b,a,matvecdata,matdata,vecdata)
end

"""
"""
function assemble_matrix(a::Assembler,cellmat,cellidsrows,cellidscols)
  A = allocate_matrix(a,cellidsrows,cellidscols)
  assemble_matrix!(A,a,cellmat,cellidsrows,cellidscols)
  A
end

"""
"""
function assemble_vector(a::Assembler,cellvec,cellids)
  b = allocate_vector(a,cellids)
  assemble_vector!(b,a,cellvec,cellids)
  b
end

"""
"""
function assemble_matrix_and_vector(a::Assembler,matvecdata,matdata,vecdata)
  A, b = allocate_matrix_and_vector(a,matvecdata,matdata,vecdata)
  assemble_matrix_and_vector!(A,b,a,matvecdata,matdata,vecdata)
  (A, b)
end

function assemble_matrix_and_vector(a::Assembler,matvecdata)
  matdata = ([],[],[])
  vecdata = ([],[])
  assemble_matrix_and_vector(a,matvecdata,matdata,vecdata)
end

function _rearange_cell_ids(matvecdata,matdata,vecdata)

  matvec1, rows1, cols1 = matvecdata
  mat2, rows2, cols2 = matdata
  vec3, rows3 = vecdata

  matrows = vcat(rows1,rows2)
  matcols = vcat(cols1,cols2)
  vecrows = vcat(rows1,rows3)

  (matrows, matcols, vecrows)
end

"""
"""
function test_assembler(a::Assembler,matvecdata,matdata,vecdata)
  trial_fesp = get_trial(a)
  test_fesp = get_test(a)
  A = allocate_matrix(a,matdata...)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  assemble_matrix!(A,a,matdata...)
  assemble_matrix_add!(A,a,matdata...)
  A = assemble_matrix(a,matdata...)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  b = allocate_vector(a,vecdata...)
  @test num_free_dofs(test_fesp) == length(b)
  assemble_vector!(b,a,vecdata...)
  assemble_vector_add!(b,a,vecdata...)
  b = assemble_vector(a,vecdata...)
  @test num_free_dofs(test_fesp) == length(b)
  A, b = allocate_matrix_and_vector(a,matvecdata,matdata,vecdata)
  assemble_matrix_and_vector!(A,b,a,matvecdata,matdata,vecdata)
  assemble_matrix_and_vector_add!(A,b,a,matvecdata,matdata,vecdata)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  @test num_free_dofs(test_fesp) == length(b)
  A, b = assemble_matrix_and_vector(a,matvecdata,matdata,vecdata)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  @test num_free_dofs(test_fesp) == length(b)
  strategy = get_assembly_strategy(a)
  @test isa(strategy,AssemblyStrategy)
end

# This is an extended interface that only make sense for assemblers that build sparse matrices
# (e.g. not for matrix free assemblers)

"""
"""
abstract type SparseMatrixAssembler <: Assembler end

"""
"""
function get_matrix_type(a::SparseMatrixAssembler)
  @abstractmethod
end

"""
"""
function get_vector_type(a::SparseMatrixAssembler)
  @abstractmethod
end

function allocate_vector(a::SparseMatrixAssembler,term_to_cellidsrows)
  n = num_free_dofs(a.test)
  allocate_vector(get_vector_type(a),n)
end

function assemble_vector!(b,a::SparseMatrixAssembler,term_to_cellvec,term_to_cellidsrows)
  fill_entries!(b,zero(eltype(b)))
  assemble_vector_add!(b,a,term_to_cellvec,term_to_cellidsrows)
end

"""
"""
function count_matrix_nnz_coo(a::SparseMatrixAssembler,term_to_cellidsrows, term_to_cellidscols)
  @abstractmethod
end

function count_matrix_nnz_coo(a::SparseMatrixAssembler,term_to_cellmat,term_to_cellidsrows, term_to_cellidscols)
  count_matrix_nnz_coo(a,term_to_cellidsrows, term_to_cellidscols)
end

"""
"""
function fill_matrix_coo_symbolic!(I,J,a::SparseMatrixAssembler,term_to_cellidsrows, term_to_cellidscols)
  @abstractmethod
end

function fill_matrix_coo_symbolic!(I,J,a::SparseMatrixAssembler,term_to_cellmat,term_to_cellidsrows, term_to_cellidscols)
  fill_matrix_coo_symbolic!(I,J,a,term_to_cellidsrows, term_to_cellidscols)
end

function allocate_matrix(a::SparseMatrixAssembler, term_to_cellidsrows, term_to_cellidscols)
  n = count_matrix_nnz_coo(a,term_to_cellidsrows,term_to_cellidscols)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  fill_matrix_coo_symbolic!(I,J,a,term_to_cellidsrows,term_to_cellidscols)
  m = num_free_dofs(a.test)
  n = num_free_dofs(a.trial)
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  sparse_from_coo(get_matrix_type(a),I,J,V,m,n)
end

function assemble_matrix!(
  mat,a::SparseMatrixAssembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
  z = zero(eltype(mat))
  fill_entries!(mat,z)
  assemble_matrix_add!(mat,a,term_to_cellmat,term_to_cellidsrows,term_to_cellidscols)
end

"""
"""
function fill_matrix_coo_numeric!(
  I,J,V,a::SparseMatrixAssembler,term_to_cellmat,term_to_cellidsrows, term_to_cellidscols,n=0)
  @abstractmethod
end

function assemble_matrix(
  a::SparseMatrixAssembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)

  n = count_matrix_nnz_coo(a,term_to_cellidsrows,term_to_cellidscols)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  fill_matrix_coo_numeric!(I,J,V,a,term_to_cellmat,term_to_cellidsrows,term_to_cellidscols)
  m = num_free_dofs(a.test)
  n = num_free_dofs(a.trial)
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  sparse_from_coo(get_matrix_type(a),I,J,V,m,n)
end

function assemble_matrix_and_vector!(A,b,a::SparseMatrixAssembler, matvecdata, matdata, vecdata)
  fill_entries!(A,zero(eltype(A)))
  fill_entries!(b,zero(eltype(b)))
  assemble_matrix_and_vector_add!(A,b,a,matvecdata, matdata, vecdata)
  A, b
end

"""
"""
function fill_matrix_and_vector_coo_numeric!(I,J,V,b,a::SparseMatrixAssembler,matvecdata, matdata, vecdata,n=0)
  @abstractmethod
end

function assemble_matrix_and_vector(a::SparseMatrixAssembler, matvecdata, matdata, vecdata)

  term_to_cellidsrows, term_to_cellidscols,  =  _rearange_cell_ids(matvecdata,matdata,vecdata)
  n = count_matrix_nnz_coo(a,term_to_cellidsrows,term_to_cellidscols)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  b = allocate_vector(a,vecdata...)

  fill_matrix_and_vector_coo_numeric!(I,J,V,b,a, matvecdata, matdata, vecdata)

  m = num_free_dofs(a.test)
  n = num_free_dofs(a.trial)
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  A = sparse_from_coo(get_matrix_type(a),I,J,V,m,n)

  A, b
end

function test_sparse_matrix_assembler(a::SparseMatrixAssembler,matvecdata,matdata,vecdata)
  test_assembler(a,matvecdata,matdata,vecdata)
  _ = get_matrix_type(a)
  _ = get_vector_type(a)
end

