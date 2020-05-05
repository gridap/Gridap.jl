
# The object returned by get_cell_basis has to implement the following trait

FECellBasisStyle(::Type{T}) where T = Val{false}()

FECellBasisStyle(cell_basis) = FECellBasisStyle(typeof(cell_basis))

"""
"""
function is_a_fe_cell_basis(cell_basis)
  v = FECellBasisStyle(cell_basis)
  get_val_parameter(v)
end

function is_a_fe_cell_basis(::Type)
  @unreachable "is_a_fe_cell_basis cannot be called on types"
end

"""
"""
abstract type FESpace <: GridapType end

# Minimal FE interface (used by FEOperator)

"""
"""
function num_free_dofs(f::FESpace)
  @abstractmethod
end

"""
"""
function FEFunction(fe::FESpace, free_values)
  @abstractmethod
end

function EvaluationFunction(fe::FESpace, free_values)
  FEFunction(fe,free_values)
end

"""
"""
function zero_free_values(fs::FESpace)
  @abstractmethod
end

"""
"""
function Base.zero(f::FESpace)
  free_values = zero_free_values(f)
  FEFunction(f,free_values)
end

# Extended FEInterface used by FEOperatorFromTerms and Assemblers

"""
"""
function get_cell_basis(f::FESpace)
  @abstractmethod
end

"""
"""
function constraint_style(::Type{<:FESpace})
  @abstractmethod
end

constraint_style(f::T) where T<:FESpace = constraint_style(T)

"""
"""
function get_constraint_kernel_matrix_cols(f::FESpace)
  @abstractmethod
end

"""
"""
function get_constraint_kernel_matrix_rows(f::FESpace)
  @abstractmethod
end

"""
"""
function get_constraint_kernel_vector(f::FESpace)
  @abstractmethod
end

"""
"""
function test_fe_space(f::FESpace)
  free_values = zero_free_values(f)
  fe_function = FEFunction(f,free_values)
  test_fe_function(fe_function)
  fe_basis = get_cell_basis(f)
  @test isa(has_constraints(f),Bool)
  @test isa(has_constraints(typeof(f)),Bool)
end

function test_fe_space(f::FESpace,matvecdata,matdata,vecdata)
  test_fe_space(f)

  cellmat, cellidsrows, cellidscols = matdata
  cm = apply_constraints_matrix_cols(f,cellmat,cellidscols)
  if ! has_constraints(f)
    @test cm === cellmat
  end
  cm = apply_constraints_matrix_rows(f,cellmat,cellidsrows)
  if ! has_constraints(f)
    @test cm === cellmat
  end

  cellvec, cellidsrows = vecdata
  cv = apply_constraints_vector(f,cellvec,cellidsrows)
  if ! has_constraints(f)
    @test cv === cellvec
  end

  cellmatvec, cellidsrows, cellidscols = matvecdata
  cmv = apply_constraints_matrix_and_vector_cols(f,cellmatvec,cellidscols)
  if ! has_constraints(f)
    @test cmv === cellmatvec
  end
  cmv = apply_constraints_matrix_and_vector_rows(f,cellmatvec,cellidsrows)
  if ! has_constraints(f)
    @test cmv === cellmatvec
  end

end

# API

"""
"""
function has_constraints(::Type{T}) where T <:FESpace
  v = constraint_style(T)
  get_val_parameter(v)
end

has_constraints(f::T) where T<:FESpace = has_constraints(T)

"""
"""
function apply_constraints_matrix_cols(f::FESpace,cellmat,cellids)
  _apply_constraints_matrix_cols(constraint_style(f),f,cellmat,cellids)
end

function _apply_constraints_matrix_cols(::Val{false},f,cellmat,cellids)
  cellmat
end

function _apply_constraints_matrix_cols(::Val{true},f,cellmat,cellids)
  k = get_constraint_kernel_matrix_cols(f)
  apply(k,cellmat,cellids)
end

"""
"""
function apply_constraints_matrix_rows(f::FESpace,cellmat,cellids)
  _apply_constraints_matrix_rows(constraint_style(f),f,cellmat,cellids)
end

function _apply_constraints_matrix_rows(::Val{false},f,cellmat,cellids)
  cellmat
end

function _apply_constraints_matrix_rows(::Val{true},f,cellmat,cellids)
  k = get_constraint_kernel_matrix_rows(f)
  apply(k,cellmat,cellids)
end

"""
"""
function apply_constraints_vector(f::FESpace,cellvec,cellids)
  _apply_constraints_vector(constraint_style(f),f,cellvec,cellids)
end

function _apply_constraints_vector(::Val{false},f,cellvec,cellids)
  cellvec
end

function _apply_constraints_vector(::Val{true},f,cellvec,cellids)
  k = get_constraint_kernel_vector(f)
  apply(k,cellvec,cellids)
end

"""
"""
function apply_constraints_matrix_and_vector_cols(f::FESpace,cellmatvec,cellids)
  _apply_constraints_matrix_and_vector_cols(constraint_style(f),f,cellmatvec,cellids)
end

function _apply_constraints_matrix_and_vector_cols(::Val{false},f,cellmatvec,cellids)
  cellmatvec
end

function _apply_constraints_matrix_and_vector_cols(::Val{true},f,cellmatvec,cellids)
  kmat = get_constraint_kernel_matrix_cols(f)
  k = MatKernel(kmat)
  apply(k,cellmatvec,cellids)
end

"""
"""
function apply_constraints_matrix_and_vector_rows(f::FESpace,cellmatvec,cellids)
  _apply_constraints_matrix_and_vector_rows(constraint_style(f),f,cellmatvec,cellids)
end

function _apply_constraints_matrix_and_vector_rows(::Val{false},f,cellmatvec,cellids)
  cellmatvec
end

function _apply_constraints_matrix_and_vector_rows(::Val{true},f,cellmatvec,cellids)
  kmat = get_constraint_kernel_matrix_rows(f)
  kvec = get_constraint_kernel_vector(f)
  k = MatVecKernel(kmat,kvec)
  apply(k,cellmatvec,cellids)
end

# Helpers

struct MatVecKernel{A<:Kernel,B<:Kernel} <: Kernel
  kmat::A
  kvec::B
end

function kernel_cache(k::MatVecKernel,matvec,cellid)
  mat, vec = matvec
  cmat = kernel_cache(k.kmat,mat,cellid)
  cvec = kernel_cache(k.kvec,vec,cellid)
  (cmat, cvec)
end

function kernel_return_type(k::MatVecKernel,matvec,cellid)
  mat, vec = matvec
  A = kernel_return_type(k.kmat,mat,cellid)
  B = kernel_return_type(k.kvec,vec,cellid)
  Tuple{A,B}
end

@inline function apply_kernel!(cache,k::MatVecKernel,matvec,cellid)
  cmat, cvec = cache
  mat, vec = matvec
  a = apply_kernel!(cmat,k.kmat,mat,cellid)
  b = apply_kernel!(cvec,k.kvec,vec,cellid)
  (a,b)
end

struct MatKernel{A<:Kernel} <: Kernel
  kmat::A
end

function kernel_cache(k::MatKernel,matvec,cellid)
  mat, vec = matvec
  cmat = kernel_cache(k.kmat,mat,cellid)
  cmat
end

function kernel_return_type(k::MatKernel,matvec,cellid)
  mat, vec = matvec
  A = kernel_return_type(k.kmat,mat,cellid)
  B = typeof(vec)
  Tuple{A,B}
end

@inline function apply_kernel!(cmat,k::MatKernel,matvec,cellid)
  mat, vec = matvec
  a = apply_kernel!(cmat,k.kmat,mat,cellid)
  (a,vec)
end
