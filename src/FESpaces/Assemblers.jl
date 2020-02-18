
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
function assemble_vector!(b,a::Assembler,cellvec,cellids)
  @abstractmethod
end

"""
"""
function assemble_matrix_and_vector!(A,b,a::Assembler,matvecdata,matdata,vecdata)
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
  A = assemble_matrix(a,matdata...)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  b = allocate_vector(a,vecdata...)
  @test num_free_dofs(test_fesp) == length(b)
  assemble_vector!(b,a,vecdata...)
  b = assemble_vector(a,vecdata...)
  @test num_free_dofs(test_fesp) == length(b)
  A, b = allocate_matrix_and_vector(a,matvecdata,matdata,vecdata)
  assemble_matrix_and_vector!(A,b,a,matvecdata,matdata,vecdata)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  @test num_free_dofs(test_fesp) == length(b)
  A, b = assemble_matrix_and_vector(a,matvecdata,matdata,vecdata)
  @test num_free_dofs(trial_fesp) == size(A,2)
  @test num_free_dofs(test_fesp) == size(A,1)
  @test num_free_dofs(test_fesp) == length(b)
end

