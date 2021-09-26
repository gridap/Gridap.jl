
function FESpaces._change_argument(
  op::typeof(jacobian),f,trian,uh::MultiFieldFEFunction)

  U = get_fe_space(uh)
  function g(cell_u)
    single_fields = GenericCellField[]
    nfields = length(U.spaces)
    for i in 1:nfields
      cell_values_field = lazy_map(a->a.array[i],cell_u)
      cf = CellField(U.spaces[i],cell_values_field)
      cell_data = lazy_map(BlockMap((1,nfields),i),get_data(cf))
      uhi = GenericCellField(cell_data,get_triangulation(cf),DomainStyle(cf))
      push!(single_fields,uhi)
    end
    xh = MultiFieldCellField(single_fields)
    cell_grad = f(xh)
    get_contribution(cell_grad,trian)
  end
  g
end

function FESpaces._change_argument(
  op::typeof(gradient),f,trian,uh::MultiFieldFEFunction)

  U = get_fe_space(uh)
  function g(cell_u)
    single_fields = GenericCellField[]
    nfields = length(U.spaces)
    for i in 1:nfields
      cell_values_field = lazy_map(a->a.array[i],cell_u)
      cf = CellField(U.spaces[i],cell_values_field)
      cell_data = lazy_map(BlockMap(nfields,i),get_data(cf))
      uhi = GenericCellField(cell_data,get_triangulation(cf),DomainStyle(cf))
      push!(single_fields,uhi)
    end
    xh = MultiFieldCellField(single_fields)
    cell_grad = f(xh)
    get_contribution(cell_grad,trian)
  end
  g
end

function Algebra.hessian(f::Function,uh::MultiFieldFEFunction)
  @notimplemented
end

#function FESpaces._change_argument(
#  op::typeof(hessian),f,trian,uh::MultiFieldFEFunction)
#
#  U = get_fe_space(uh)
#  function g(cell_u)
#    single_fields = GenericCellField[]
#    nfields = length(U.spaces)
#    for i in 1:nfields
#      cell_values_field = lazy_map(a->a.array[i],cell_u)
#      cf = CellField(U.spaces[i],cell_values_field)
#      cell_data = lazy_map(BlockMap((nfields,nfields),i),get_data(cf))
#      uhi = GenericCellField(cell_data,get_triangulation(cf),DomainStyle(cf))
#      push!(single_fields,uhi)
#    end
#    xh = MultiFieldCellField(single_fields)
#    cell_grad = f(xh)
#    get_contribution(cell_grad,trian)
#  end
#  g
#end
