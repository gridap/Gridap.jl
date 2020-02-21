
"""
"""
struct SingleFieldFEFunction <: CellField
  array
  cell_vals
  free_values
  dirichlet_values
  fe_space
  @doc """
  """
  function SingleFieldFEFunction(
    array::AbstractArray{<:Field},
    cell_vals::AbstractArray{<:AbstractArray},
    free_values::AbstractVector,
    dirichlet_values::AbstractVector,
    fe_space::SingleFieldFESpace)
    new(array,cell_vals,free_values,dirichlet_values,fe_space)
  end
end

FEFunctionStyle(::Type{<:SingleFieldFEFunction}) = Val{true}()

get_array(f::SingleFieldFEFunction) = f.array

get_free_values(f::SingleFieldFEFunction) = f.free_values

get_dirichlet_values(f::SingleFieldFEFunction) = f.dirichlet_values

get_fe_space(f::SingleFieldFEFunction) = f.fe_space

get_cell_map(f::SingleFieldFEFunction) = get_cell_map(f.fe_space)

get_cell_values(f::SingleFieldFEFunction) = f.cell_vals

