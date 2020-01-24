
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
function allocate_matrix(a::Assembler,cellidsrows,cellidscols)
  @abstractmethod
end

"""
"""
function assemble_matrix!(A,a::Assembler,cellmat,cellidsrows,cellidscols)
  @abstractmethod
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
function allocate_vector(a::Assembler,cellidsrows)
  @abstractmethod
end

"""
"""
function assemble_vector!(b,a::Assembler,cellvec,cellids)
  @abstractmethod
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
function allocate_matrix_and_vector(a::Assembler,cellidsrows,cellidscols)
  A = allocate_matrix(a,cellidsrows,cellidscols)
  b = allocate_vector(a,cellidsrows)
  (A,b)
end

"""
"""
function assemble_matrix_and_vector!(A,b,a::Assembler,cellmat,cellvec,cellidsrows,cellidscols)
  assemble_matrix!(A,a,cellmat,cellidsrows,cellidscols)
  assemble_vector!(b,a,cellvec,cellidsrows)
  (A,b)
end

"""
"""
function assemble_matrix_and_vector(a::Assembler,cellmat,cellvec,cellidsrows,cellidscols)
  A, b = allocate_matrix_and_vector(a,cellidsrows,cellidscols)
  assemble_matrix_and_vector!(A,b,cellmat,cellvec,cellidsrows,cellidscols)
  (A,b)
end

"""
"""
function test_assembler(a::Assembler,cellmat,cellvec,cellidsrows,cellidscols)
  trial_fesp = get_trial(a)
  test_fesp = get_test(a)
  A = allocate_matrix(a,cellidsrows,cellidscols)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  assemble_matrix!(A,a,cellmat,cellidsrows,cellidscols)
  b = allocate_vector(a,cellidsrows)
  @test num_free_dofs(test_fesp) == length(b)
  assemble_vector!(b,a,cellvec,cellidsrows)
end

