
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
  function fun(uh)
    trian_i = get_triangulation(uh)
    if have_compatible_domains(trian_i,trian) ||
      have_compatible_domains(trian_i,get_background_triangulation(trian)) ||
      Geometry.is_included(trian,trian_i)
      cell_dofs = get_cell_dof_values(uh,trian)
    else
      cell_dofs_i = get_cell_dof_values(uh)
      T = eltype(eltype(cell_dofs_i))
      cell_dofs = Fill(T[],num_cells(trian))
    end
    cell_dofs
  end
  uhs = f.single_fe_functions
  blocks = [ fun(uh) for uh in uhs ]
  bids = [ (i,) for i in 1:length(uhs)]
  bsize = (length(uhs),)
  cell_axes_blocks = map(i->lazy_map(axes,i),blocks)
  cell_axes = lazy_map(_multifield_axes_dofs,cell_axes_blocks...)
  lazy_map(BlockArrayCooMap(bsize,bids),cell_axes,blocks...)
end

function FESpaces.get_cell_dof_values(f::MultiFieldFEFunction,trian::SkeletonTriangulation)
  cell_values_plus = get_cell_dof_values(f,trian.plus)
  cell_values_minus = get_cell_dof_values(f,trian.minus)
  cell_axes_plus = lazy_map(axes,cell_values_plus)
  cell_axes_minus = lazy_map(axes,cell_values_minus)
  cell_axes = lazy_map(cell_axes_plus,cell_axes_minus) do axp, axm
    (append_ranges([axp[1],axm[1]]),)
  end
  lazy_map(BlockArrayCooMap((2,),[(1,),(2,)]),cell_axes,cell_values_plus,cell_values_minus)
end

"""
    num_fields(m::MultiFieldFEFunction)
"""
num_fields(m::MultiFieldFEFunction) = length(m.single_fe_functions)
Base.iterate(m::MultiFieldFEFunction) = iterate(m.single_fe_functions)
Base.iterate(m::MultiFieldFEFunction,state) = iterate(m.single_fe_functions,state)
Base.getindex(m::MultiFieldFEFunction,field_id::Integer) = m.single_fe_functions[field_id]
