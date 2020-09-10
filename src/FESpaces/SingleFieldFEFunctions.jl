
"""
"""
struct SingleFieldFEFunction{T<:CellField} <: CellField
  cell_field::T
  cell_vals
  free_values
  dirichlet_values
  fe_space
  @doc """
  """
  function SingleFieldFEFunction(
    cell_field::CellField,
    cell_vals::AbstractArray{<:AbstractArray},
    free_values::AbstractVector,
    dirichlet_values::AbstractVector,
    fe_space::SingleFieldFESpace)

    T = typeof(cell_field)
    new{T}(cell_field,cell_vals,free_values,dirichlet_values,fe_space)
  end
end

FEFunctionStyle(::Type{<:SingleFieldFEFunction}) = Val{true}()

get_array(f::SingleFieldFEFunction) = get_array(f.cell_field)

get_free_values(f::SingleFieldFEFunction) = f.free_values

get_dirichlet_values(f::SingleFieldFEFunction) = f.dirichlet_values

get_fe_space(f::SingleFieldFEFunction) = f.fe_space

get_cell_map(f::SingleFieldFEFunction) = get_cell_map(f.cell_field)

get_cell_values(f::SingleFieldFEFunction) = f.cell_vals

CellData.MetaSizeStyle(::Type{SingleFieldFEFunction{T}}) where T = MetaSizeStyle(T)

CellData.get_cell_axes(a::SingleFieldFEFunction) = get_cell_axes(a.cell_field)

CellData.get_memo(a::SingleFieldFEFunction) = get_memo(a.cell_field)

#TODO this is a bit hacky
CellData._compose_cell_fields(f::SingleFieldFEFunction,ϕ) = f.cell_field∘ϕ

Fields.gradient(f::SingleFieldFEFunction) = gradient(f.cell_field)
