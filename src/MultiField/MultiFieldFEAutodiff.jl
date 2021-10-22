
# TO IMPROVE .... Perhaps via a Map?
# Do we have already this function in Gridap?
function _get_cell_dofs_field_offsets(uh::MultiFieldFEFunction)
  U = get_fe_space(uh)
  uh_dofs = get_cell_dof_values(uh)[1]
  nfields = length(U.spaces)
  dofs_field_offsets=Vector{Int}(undef,nfields+1)
  dofs_field_offsets[1]=1
  for i in 1:nfields
    dofs_field_offsets[i+1]=dofs_field_offsets[i]+length(uh_dofs.array[i])
  end
  dofs_field_offsets
end

function FESpaces._gradient(f,uh::MultiFieldFEFunction,fuh::DomainContribution)
  terms = DomainContribution()
  U = get_fe_space(uh)
  for trian in get_domains(fuh)
    g = FESpaces._change_argument(gradient,f,trian,uh)
    cell_u = lazy_map(DensifyInnerMostBlockLevelMap(),get_cell_dof_values(uh))
    glue = get_glue(trian,Val(Geometry.num_cell_dims(get_triangulation(uh))))
    cell_id = FESpaces._compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_gradient(g,cell_u,cell_id)
    monolithic_result=cell_grad
    block_result = [] # TO-DO type unstable
    nfields = length(U.spaces)
    cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
    for i in 1:nfields
      view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
      block=lazy_map(x->view(x,view_range),monolithic_result)
      append!(block_result,[block])
    end
    cell_grad=lazy_map(BlockMap(nfields,collect(1:nfields)),block_result...)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

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
    cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
    for i in 1:nfields
      view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
      cell_values_field = lazy_map(a->view(a,view_range),cell_u)
      cf = CellField(U.spaces[i],cell_values_field)
      push!(single_fields,cf)
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
