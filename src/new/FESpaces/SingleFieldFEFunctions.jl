
"""
"""
struct SingleFieldFEFunction{A,B,C} <: CellField
  array::A
  free_values::B
  dirichlet_values::B
  fe_space::C
  @doc """
  """
  function SingleFieldFEFunction(
    array::AbstractArray{<:Field},
    free_values::AbstractVector,
    dirichlet_values::AbstractVector,
    fe_space::SingleFieldFESpace)

    A = typeof(array)
    B = typeof(free_values)
    C = typeof(fe_space)
    new{A,B,C}(array,free_values,dirichlet_values,fe_space)
  end
end

FEFunctionStyle(::Type{<:SingleFieldFEFunction}) = Val{true}()

get_array(f::SingleFieldFEFunction) = f.array

get_free_values(f::SingleFieldFEFunction) = f.free_values

get_dirichlet_values(f::SingleFieldFEFunction) = f.dirichlet_values

get_fe_space(f::SingleFieldFEFunction) = f.fe_space

get_cell_map(f::SingleFieldFEFunction) = get_cell_map(f.fe_space)

