
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

"""
"""
function num_free_dofs(f::FESpace)
  @abstractmethod
end

"""
"""
function get_cell_basis(f::FESpace)
  @abstractmethod
end

"""
"""
function FEFunction(fe::FESpace, free_values)
  @abstractmethod
end

"""
"""
function zero_free_values(::Type{T},fs::FESpace) where T
  @abstractmethod
end

function zero_free_values(fs::FESpace)
  zero_free_values(Float64,fs)
end

"""
"""
function Base.zero(f::FESpace)
  free_values = zero_free_values(f)
  FEFunction(f,free_values)
end

"""
"""
function apply_constraints_matrix_cols(f::FESpace,cellmat,cellids)
  @abstractmethod
end

"""
"""
function apply_constraints_matrix_rows(f::FESpace,cellmat,cellids)
  @abstractmethod
end

"""
"""
function apply_constraints_vector(f::FESpace,cellvec,cellids)
  @abstractmethod
end

"""
"""
function apply_constraints_matrix_and_vector_rows(f::FESpace,cellmat,cellvec,cellids)
  cm = apply_constraints_matrix_rows(f,cellmat,cellids)
  cv = apply_constraints_vector(f,cellvec,cellids)
  (cm, cv)
end

"""
"""
function test_fe_space(f::FESpace,cellmat,cellvec,cellidsrows,cellidscols)
  free_values = zero_free_values(f)
  @test eltype(zero_free_values(Int,f)) == Int
  fe_function = FEFunction(f,free_values)
  test_fe_function(fe_function)
  fe_basis = get_cell_basis(f)
  _ = apply_constraints_matrix_cols(f,cellmat,cellidscols)
  _ = apply_constraints_matrix_rows(f,cellmat,cellidsrows)
  _ = apply_constraints_vector(f,cellvec,cellidsrows)
end

