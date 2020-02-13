
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
function allocate_matrix_and_vector(a::Assembler,matvecdata,matdata,vecdata)
  (_, matrows, matcols), (_, vecrows) = _rearange_matvec_data(matvecdata,matdata,vecdata)
  A = allocate_matrix(a,matrows,matcols)
  b = allocate_vector(a,vecrows)
  (A,b)
end

"""
"""
function assemble_matrix_and_vector!(A,b,a::Assembler,matvecdata,matdata,vecdata)
  _matdata, _vecdata = _rearange_matvec_data(matvecdata,matdata,vecdata)
  assemble_matrix!(A,a,_matdata...)
  assemble_vector!(b,a,_vecdata...)
  (A,b)
end

"""
"""
function assemble_matrix_and_vector(a::Assembler,matvecdata,matdata,vecdata)
  _matdata, _vecdata = _rearange_matvec_data(matvecdata,matdata,vecdata)
  A = assemble_matrix(a,_matdata...)
  b = assemble_vector(a,_vecdata...)
  (A,b)
end

function _rearange_matvec_data(matvecdata,matdata,vecdata)

  mat1, vec1, rows1, cols1 = matvecdata
  mat2, rows2, cols2 = matdata
  vec3, rows3 = vecdata

  mat = vcat(mat1,mat2)
  matrows = vcat(rows1,rows2)
  matcols = vcat(cols1,cols2)

  vec = vcat(vec1,vec3)
  vecrows = vcat(rows1,rows3)

  ((mat,matrows,matcols), (vec,vecrows))
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

