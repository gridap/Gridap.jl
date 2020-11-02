
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

