
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
function allocate_matrix(a::Assembler)
  @abstractmethod
end

"""
"""
function assemble_matrix!(A,a::Assembler,cellmat)
  @abstractmethod
end

"""
"""
function assemble_matrix(a::Assembler,cellmat)
  A = allocate_matrix(a)
  assemble_matrix!(A,a,cellmat)
  A
end

"""
"""
function allocate_vector(a::Assembler)
  @abstractmethod
end

"""
"""
function assemble_vector!(b,a::Assembler,cellvec)
  @abstractmethod
end

"""
"""
function assemble_vector(a::Assembler,cellvec)
  b = allocate_vector(a)
  assemble_vector!(b,a,cellvec)
  b
end

"""
"""
function allocate_matrix_and_vector(a::Assembler,cellmat,cellvec)
  A = allocate_matrix(a)
  b = allocate_vector(a)
  (A,b)
end

"""
"""
function assemble_matrix_and_vector!(A,b,a::Assembler,cellmat,cellvec)
  assemble_matrix!(A,a,cellmat)
  assemble_vector!(b,a,cellvec)
  (A,b)
end

"""
"""
function assemble_matrix_and_vector(a::Assembler,cellmat,cellvec)
  A, b = allocate_matrix_and_vector(a)
  assemble_matrix_and_vector!(A,b,cellmat,cellvec)
  (A,b)
end

"""
"""
function test_assembler(a::Assembler,cellmat,cellvec)
  trial_fesp = get_trial(a)
  test_fesp = get_test(a)
  A = allocate_matrix(a)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  assemble_matrix!(A,a,cellmat)
  b = allocate_vector(a)
  @test num_free_dofs(test_fesp) == length(b)
  assemble_matrix!(b,a,cellvec)
end

