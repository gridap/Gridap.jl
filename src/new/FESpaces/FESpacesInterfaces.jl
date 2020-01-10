
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
function get_cell_fe_basis(f::FESpace)
  @abstractmethod
end

"""
"""
function FEFunction(fe::FESpace, free_values)
  @abstractmethod
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

"""
"""
function apply_constraints_matrix_cols(f::FESpace,cellmat)
  @abstractmethod
end

"""
"""
function apply_constraints_matrix_rows(f::FESpace,cellmat)
  @abstractmethod
end

"""
"""
function apply_constraints_vector(f::FESpace,cellvec)
  @abstractmethod
end

"""
"""
function apply_constraints_matrix_and_vector_rows(f::FESpace,cellmat,cellvec)
  cm = apply_constraints_matrix_rows(f,cellmat)
  cv = apply_constraints_vector(f,cellvec)
  (cm, cv)
end

"""
"""
function test_fe_space(f::FESpace,cellmat,cellvec)
  free_values = zero_free_values(f)
  fe_function = FEFunction(f,free_values)
  fe_basis = get_cell_fe_basis(f)
  _ = apply_constraints_matrix_cols(f,cellmat)
  _ = apply_constraints_matrix_rows(f,cellmat)
  _ = apply_constraints_vector(f,cellmat)
end

