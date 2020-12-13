
"""
    struct MultiFieldFEFunction <: CellField
      # private fields
    end
"""
struct MultiFieldFEFunction{T<:MultiFieldCellField} <: FEFunction
  single_fe_functions::Vector{<:SingleFieldFEFunction}
  free_values::AbstractArray
  cell_dof_values::AbstractArray
  fe_space::MultiFieldFESpace
  multi_cell_field::T

  function MultiFieldFEFunction(
    free_values::AbstractVector,
    space::MultiFieldFESpace,
    single_fe_functions::Vector{<:SingleFieldFEFunction})

    multi_cell_field = MultiFieldCellField(map(i->i.cell_field,single_fe_functions))
    T = typeof(multi_cell_field)

    blocks = map(get_cell_dof_values,single_fe_functions)
    cell_axes_blocks = map(i->lazy_map(axes,i),blocks)
    cell_axes = lazy_map(_multifield_axes_dofs,cell_axes_blocks...)
    blockids = [(i,) for i in 1:length(single_fe_functions)]
    bsize = (length(single_fe_functions),)
    cell_dof_values = lazy_map(BlockArrayCooMap(bsize,blockids),cell_axes,blocks...)

    new{T}(
      single_fe_functions,
      free_values,
      cell_dof_values,
      space,
      multi_cell_field)
  end
end

CellData.get_data(f::MultiFieldFEFunction) = get_data(f.multi_cell_field)
CellData.get_triangulation(f::MultiFieldFEFunction) = get_triangulation(f.multi_cell_field)
CellData.DomainStyle(::Type{MultiFieldFEFunction{T}}) where T = DomainStyle(T)
FESpaces.get_free_values(f::MultiFieldFEFunction) = f.free_values
FESpaces.get_fe_space(f::MultiFieldFEFunction) = f.fe_space
FESpaces.get_cell_dof_values(f::MultiFieldFEFunction) = f.cell_dof_values

"""
    num_fields(m::MultiFieldFEFunction)
"""
num_fields(m::MultiFieldFEFunction) = length(m.single_fe_functions)
Base.iterate(m::MultiFieldFEFunction) = iterate(m.single_fe_functions)
Base.iterate(m::MultiFieldFEFunction,state) = iterate(m.single_fe_functions,state)
Base.getindex(m::MultiFieldFEFunction,field_id::Integer) = m.single_fe_functions[field_id]
