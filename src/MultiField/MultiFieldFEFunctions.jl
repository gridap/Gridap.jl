
"""
    struct MultiFieldFEFunction <: CellField
      # private fields
    end
"""
struct MultiFieldFEFunction{T<:MultiFieldCellField} <: FEFunction
  single_fe_functions::Vector{<:SingleFieldFEFunction}
  free_values::AbstractArray
  fe_space::MultiFieldFESpace
  multi_cell_field::T

  function MultiFieldFEFunction(
    free_values::AbstractVector,
    space::MultiFieldFESpace,
    single_fe_functions::Vector{<:SingleFieldFEFunction})

    multi_cell_field = MultiFieldCellField(map(i->i.cell_field,single_fe_functions))
    T = typeof(multi_cell_field)

    new{T}(
      single_fe_functions,
      free_values,
      space,
      multi_cell_field)
  end
end

CellData.get_data(f::MultiFieldFEFunction) = get_data(f.multi_cell_field)
CellData.get_triangulation(f::MultiFieldFEFunction) = get_triangulation(f.multi_cell_field)
CellData.DomainStyle(::Type{MultiFieldFEFunction{T}}) where T = DomainStyle(T)
FESpaces.get_free_dof_values(f::MultiFieldFEFunction) = f.free_values
FESpaces.get_fe_space(f::MultiFieldFEFunction) = f.fe_space

function FESpaces.get_cell_dof_values(f::MultiFieldFEFunction)
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get the DOF values
  on top of the corresponding cells.
  """
  trians = map(get_triangulation,f.fe_space.spaces)
  trian = first(trians)
  @check all(map(t->have_compatible_domains(t,trian),trians)) msg
  get_cell_dof_values(f,trian)
end

function FESpaces.get_cell_dof_values(f::MultiFieldFEFunction,trian::Triangulation)
  uhs = f.single_fe_functions
  blockmask = [ _can_be_restricted_to(get_triangulation(uh),trian) for uh in uhs ]
  active_block_ids = findall(blockmask)
  active_block_data = Any[ get_cell_dof_values(uhs[i],trian) for i in active_block_ids ]
  nblocks = length(uhs)
  lazy_map(BlockMap(nblocks,active_block_ids),active_block_data...)
end

function FESpaces.get_cell_dof_values(f::MultiFieldFEFunction,trian::SkeletonTriangulation)
  cell_values_plus = get_cell_dof_values(f,trian.plus)
  cell_values_minus = get_cell_dof_values(f,trian.minus)
  lazy_map(BlockMap(2,[1,2]),cell_values_plus,cell_values_minus)
end

"""
    num_fields(m::MultiFieldFEFunction)
"""
num_fields(m::MultiFieldFEFunction) = length(m.single_fe_functions)
Base.iterate(m::MultiFieldFEFunction) = iterate(m.single_fe_functions)
Base.iterate(m::MultiFieldFEFunction,state) = iterate(m.single_fe_functions,state)
Base.getindex(m::MultiFieldFEFunction,field_id::Integer) = m.single_fe_functions[field_id]
