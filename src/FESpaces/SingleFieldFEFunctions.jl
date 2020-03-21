
"""
"""
struct SingleFieldFEFunction{R} <: CellField
  array
  cell_vals
  free_values
  dirichlet_values
  fe_space
  ref_style::Val{R}
  @doc """
  """
  function SingleFieldFEFunction(
    array::AbstractArray{<:Field},
    cell_vals::AbstractArray{<:AbstractArray},
    free_values::AbstractVector,
    dirichlet_values::AbstractVector,
    fe_space::SingleFieldFESpace)

    ref_style = RefStyle(get_cell_dof_basis(fe_space))
    R = get_val_parameter(ref_style)
    new{R}(array,cell_vals,free_values,dirichlet_values,fe_space,ref_style)
  end
end

FEFunctionStyle(::Type{<:SingleFieldFEFunction}) = Val{true}()

get_array(f::SingleFieldFEFunction) = f.array

get_free_values(f::SingleFieldFEFunction) = f.free_values

get_dirichlet_values(f::SingleFieldFEFunction) = f.dirichlet_values

get_fe_space(f::SingleFieldFEFunction) = f.fe_space

get_cell_map(f::SingleFieldFEFunction) = get_cell_map(f.fe_space)

get_cell_values(f::SingleFieldFEFunction) = f.cell_vals

RefStyle(::Type{SingleFieldFEFunction{R}}) where R = Val{R}()
