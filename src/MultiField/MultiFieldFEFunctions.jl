
"""
    struct MultiFieldFEFunction <: CellField
      # private fields
    end
"""
struct MultiFieldFEFunction{T<:MultiFieldCellField} <: CellField
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


    f(i) = get_cell_axes(get_fe_space(i))
    cell_axes = create_array_of_blocked_axes(map(f,single_fe_functions)...)
    blocks = Tuple(map(get_cell_values,single_fe_functions))
    blockids = [(i,) for i in 1:length(single_fe_functions)]
    cell_dof_values = VectorOfBlockArrayCoo(blocks,blockids,cell_axes)

    new{T}(
      single_fe_functions,
      free_values,
      cell_dof_values,
      space,
      multi_cell_field)
  end
end

#function MultiFieldFEFunction(
#  space::MultiFieldFESpace,
#  blocks::Vector{<:SingleFieldFEFunction})
#  fv = zero_free_values(space)
#  xh0 = MultiFieldFEFunction(fv,space,blocks)
#end

Arrays.get_array(a::MultiFieldFEFunction) = get_array(a.multi_cell_field)

CellData.get_memo(a::MultiFieldFEFunction) = get_memo(a.multi_cell_field)

CellData.get_cell_map(a::MultiFieldFEFunction) = get_cell_map(a.multi_cell_field)

CellData.get_cell_axes(a::MultiFieldFEFunction) = get_cell_axes(a.multi_cell_field)

CellData.RefStyle(::Type{<:MultiFieldFEFunction{T}}) where T = RefStyle(T)

CellData.MetaSizeStyle(::Type{<:MultiFieldFEFunction{T}}) where T = MetaSizeStyle(T)

FESpaces.FEFunctionStyle(::Type{<:MultiFieldFEFunction}) = Val{true}()

FESpaces.get_free_values(f::MultiFieldFEFunction) = f.free_values

FESpaces.get_fe_space(f::MultiFieldFEFunction) = f.fe_space

"""
    num_fields(m::MultiFieldFEFunction)
"""
num_fields(m::MultiFieldFEFunction) = length(m.single_fe_functions)

Base.iterate(m::MultiFieldFEFunction) = iterate(m.single_fe_functions)

Base.iterate(m::MultiFieldFEFunction,state) = iterate(m.single_fe_functions,state)

Base.getindex(m::MultiFieldFEFunction,field_id::Integer) = m.single_fe_functions[field_id]

function Geometry.restrict(a::MultiFieldFEFunction,trian::Triangulation)
  f = (ai) -> restrict(ai,trian)
  blocks = map(f,a.single_fe_functions)
  blocks
end

function FESpaces.get_cell_values(f::MultiFieldFEFunction)
  f.cell_dof_values
end

