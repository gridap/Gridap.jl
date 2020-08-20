
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
function allocate_matrix(a::Assembler,matdata)
  @abstractmethod
end

"""
"""
function allocate_vector(a::Assembler,vecdata)
  @abstractmethod
end

"""
"""
function allocate_matrix_and_vector(a::Assembler,data)
  @abstractmethod
end

"""
"""
function assemble_matrix!(A,a::Assembler,matdata)
  @abstractmethod
end

"""
"""
function assemble_matrix_add!(mat,a::Assembler, matdata)
  @abstractmethod
end

"""
"""
function assemble_vector!(b,a::Assembler,vecdata)
  @abstractmethod
end

"""
"""
function assemble_vector_add!(b,a::Assembler,vecdata)
  @abstractmethod
end

"""
"""
function assemble_matrix_and_vector!(A,b,a::Assembler, data)
  @abstractmethod
end

function assemble_matrix_and_vector_add!(A,b,a::Assembler, data)
  @abstractmethod
end

"""
"""
function assemble_matrix(a::Assembler,matdata)
  A = allocate_matrix(a,matdata)
  assemble_matrix!(A,a,matdata)
  A
end

"""
"""
function assemble_vector(a::Assembler,vecdata)
  b = allocate_vector(a,vecdata)
  assemble_vector!(b,a,vecdata)
  b
end

"""
"""
function assemble_matrix_and_vector(a::Assembler,data)
  A, b = allocate_matrix_and_vector(a,data)
  assemble_matrix_and_vector!(A,b,a,data)
  (A, b)
end

"""
"""
function test_assembler(a::Assembler,matdata,vecdata,data)
  trial_fesp = get_trial(a)
  test_fesp = get_test(a)
  A = allocate_matrix(a,matdata)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  assemble_matrix!(A,a,matdata)
  assemble_matrix_add!(A,a,matdata)
  A = assemble_matrix(a,matdata)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  b = allocate_vector(a,vecdata)
  @test num_free_dofs(test_fesp) == length(b)
  assemble_vector!(b,a,vecdata)
  assemble_vector_add!(b,a,vecdata)
  b = assemble_vector(a,vecdata)
  @test num_free_dofs(test_fesp) == length(b)
  A, b = allocate_matrix_and_vector(a,data)
  assemble_matrix_and_vector!(A,b,a,data)
  assemble_matrix_and_vector_add!(A,b,a,data)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  @test num_free_dofs(test_fesp) == length(b)
  A, b = assemble_matrix_and_vector(a,data)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  @test num_free_dofs(test_fesp) == length(b)
end

# This is an extended interface that only makes sense for assemblers that build (sequential) sparse matrices
# (e.g. not for matrix free assemblers or for distributed assemblers)

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

function allocate_vector(a::SparseMatrixAssembler,vecdata)
  n = num_free_dofs(get_test(a))
  allocate_vector(get_vector_type(a),n)
end

function assemble_vector!(b,a::SparseMatrixAssembler,vecdata)
  fill_entries!(b,zero(eltype(b)))
  assemble_vector_add!(b,a,vecdata)
end

"""
"""
function count_matrix_nnz_coo(a::SparseMatrixAssembler,matdata)
  @abstractmethod
end

"""
"""
function count_matrix_and_vector_nnz_coo(a::SparseMatrixAssembler,data)
  @abstractmethod
end

"""
"""
function fill_matrix_coo_symbolic!(I,J,a::SparseMatrixAssembler,matdata,n=0)
  @abstractmethod
end

function fill_matrix_and_vector_coo_symbolic!(I,J,a::SparseMatrixAssembler,data,n=0)
  @abstractmethod
end

function allocate_matrix(a::SparseMatrixAssembler,matdata)
  n = count_matrix_nnz_coo(a,matdata)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  fill_matrix_coo_symbolic!(I,J,a,matdata)
  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  sparse_from_coo(get_matrix_type(a),I,J,V,m,n)
end

function assemble_matrix!(mat,a::SparseMatrixAssembler,matdata)
  z = zero(eltype(mat))
  fill_entries!(mat,z)
  assemble_matrix_add!(mat,a,matdata)
end

"""
"""
function fill_matrix_coo_numeric!(I,J,V,a::SparseMatrixAssembler,matdata,n=0)
  @abstractmethod
end

function assemble_matrix(a::SparseMatrixAssembler,matdata)

  n = count_matrix_nnz_coo(a,matdata)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)

  fill_matrix_coo_numeric!(I,J,V,a,matdata)

  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  sparse_from_coo(get_matrix_type(a),I,J,V,m,n)
end


function allocate_matrix_and_vector(a::SparseMatrixAssembler,data)

  n = count_matrix_and_vector_nnz_coo(a,data)

  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  fill_matrix_and_vector_coo_symbolic!(I,J,a,data)
  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  A = sparse_from_coo(get_matrix_type(a),I,J,V,m,n)

  b = allocate_vector(get_vector_type(a),m)

  A,b
end

function assemble_matrix_and_vector!(A,b,a::SparseMatrixAssembler, data)
  fill_entries!(A,zero(eltype(A)))
  fill_entries!(b,zero(eltype(b)))
  assemble_matrix_and_vector_add!(A,b,a,data)
  A, b
end

"""
"""
function fill_matrix_and_vector_coo_numeric!(I,J,V,b,a::SparseMatrixAssembler,data,n=0)
  @abstractmethod
end

function assemble_matrix_and_vector(a::SparseMatrixAssembler, data)

  n = count_matrix_and_vector_nnz_coo(a,data)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  n = num_free_dofs(get_test(a))
  b = allocate_vector(get_vector_type(a),n)

  fill_matrix_and_vector_coo_numeric!(I,J,V,b,a,data)

  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  A = sparse_from_coo(get_matrix_type(a),I,J,V,m,n)

  A, b
end

function test_sparse_matrix_assembler(a::SparseMatrixAssembler,matdata,vecdata,data)
  test_assembler(a,matdata,vecdata,data)
  _ = get_matrix_type(a)
  _ = get_vector_type(a)
end
