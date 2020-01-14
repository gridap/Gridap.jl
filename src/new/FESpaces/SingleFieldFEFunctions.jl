
"""
"""
struct SingleFieldFEFunction{A,B,C} <: FEFunction
  cell_field::A
  free_values::B
  dirichlet_values::B
  fe_space::C
  @doc """
  """
  function SingleFieldFEFunction(
    cell_field::AbstractArray{<:Field},
    free_values::AbstractVector,
    dirichlet_values::AbstractVector,
    fe_space::SingleFieldFESpace)

    A = typeof(cell_field)
    B = typeof(free_values)
    C = typeof(fe_space)
    new{A,B,C}(cell_field,free_values,dirichlet_values,fe_space)
  end
end

get_array(f::SingleFieldFEFunction) = f.cell_field

get_free_values(f::SingleFieldFEFunction) = f.free_values

get_dirichlet_values(f::SingleFieldFEFunction) = f.dirichlet_values

get_fe_space(f::SingleFieldFEFunction) = f.fe_space

# Unary operations

for op in (:+,:-)
  @eval begin
    function ($op)(a::SingleFieldFEFunction)
      apply_to_field_array(bcast($op),get_array(a))
    end
  end
end

for op in (:+,:-,:*)
  @eval begin
    function ($op)(a::SingleFieldFEFunction,b)
      fe_space = get_fe_space(a)
      cell_map = get_cell_map(fe_space)
      _b = _convert_to_integrable(b,cell_map)
      apply_to_field_array(bcast($op),get_array(a),get_array(_b))
    end
  end
end

