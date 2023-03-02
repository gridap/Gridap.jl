
# TO IMPROVE .... Perhaps via a Map?
# Do we have already this function in Gridap?
# What if there are no cells? I am assuming that there is at least one cell
# What if the number of dofs per field per cell is different among cells?
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

function _restructure_cell_grad!(
  cell_grad, uh::MultiFieldFEFunction, trian)

  monolithic_result=cell_grad
  blocks = [] # TO-DO type unstable. How can I infer the type of its entries?
  nfields = length(uh.fe_space.spaces)
  cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
  for i in 1:nfields
    view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
    block=lazy_map(x->view(x,view_range),monolithic_result)
    append!(blocks,[block])
  end
  cell_grad=lazy_map(BlockMap(nfields,collect(1:nfields)),blocks...)
  cell_grad
end

function FESpaces._gradient(f,uh::MultiFieldFEFunction,fuh::DomainContribution)
  terms = DomainContribution()
  U = get_fe_space(uh)
  for trian in get_domains(fuh)
    g = FESpaces._change_argument(gradient,f,trian,uh)
    cell_u = lazy_map(DensifyInnerMostBlockLevelMap(),get_cell_dof_values(uh))
    cell_id = FESpaces._compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_gradient(g,cell_u,cell_id)
    cell_grad = _restructure_cell_grad!(cell_grad, uh, trian)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function FESpaces._jacobian(f,uh::MultiFieldFEFunction,fuh::DomainContribution)
  terms = DomainContribution()
  U = get_fe_space(uh)
  for trian in get_domains(fuh)
    g = FESpaces._change_argument(jacobian,f,trian,uh)
    cell_u = lazy_map(DensifyInnerMostBlockLevelMap(),get_cell_dof_values(uh))
    cell_id = FESpaces._compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_jacobian(g,cell_u,cell_id)
    monolithic_result=cell_grad
    blocks        = [] # TO-DO type unstable. How can I infer the type of its entries?
    blocks_coords = Tuple{Int,Int}[]
    nfields = length(U.spaces)
    cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
    for j=1:nfields
      view_range_j=cell_dofs_field_offsets[j]:cell_dofs_field_offsets[j+1]-1
      for i=1:nfields
        view_range_i=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
        # TO-DO: depending on the residual being differentiated, we may end with
        #        blocks [i,j] full of zeros. I guess that it might desirable to early detect
        #        these zero blocks and use a touch[i,j]==false block in ArrayBlock.
        #        How can we detect that we have a zero block?
        block=lazy_map(x->view(x,view_range_i,view_range_j),monolithic_result)
        append!(blocks,[block])
        append!(blocks_coords,[(i,j)])
      end
    end
    cell_grad=lazy_map(BlockMap((nfields,nfields),blocks_coords),blocks...)
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
    cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
    for i in 1:nfields
      view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
      cell_values_field = lazy_map(a->a[view_range],cell_u)
      cf = CellField(U.spaces[i],cell_values_field)
      push!(single_fields,cf)
    end
    xh = MultiFieldCellField(single_fields)
    cell_grad = f(xh)
    cell_grad_cont_block=get_contribution(cell_grad,trian)
    bs = [cell_dofs_field_offsets[i+1]-cell_dofs_field_offsets[i] for i=1:nfields]
    lazy_map(DensifyInnerMostBlockLevelMap(),
             Fill(bs,length(cell_grad_cont_block)),
             cell_grad_cont_block)
  end
  g
end

function FESpaces._change_argument(
  op::Union{typeof(gradient),typeof(hessian)},f,trian,uh::MultiFieldFEFunction)

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

function FESpaces._hessian(f,uh::MultiFieldFEFunction,fuh::DomainContribution)
  terms = DomainContribution()
  U = get_fe_space(uh)
  for trian in get_domains(fuh)
    g = FESpaces._change_argument(hessian,f,trian,uh)
    cell_u = lazy_map(DensifyInnerMostBlockLevelMap(),get_cell_dof_values(uh))
    cell_id = FESpaces._compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_hessian(g,cell_u,cell_id)
    monolithic_result=cell_grad
    blocks        = [] # TO-DO type unstable. How can I infer the type of its entries?
    blocks_coords = Tuple{Int,Int}[]
    nfields = length(U.spaces)
    cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
    for j=1:nfields
      view_range_j=cell_dofs_field_offsets[j]:cell_dofs_field_offsets[j+1]-1
      for i=1:nfields
        view_range_i=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
        # TO-DO: depending on the residual being differentiated, we may end with
        #        blocks [i,j] full of zeros. I guess that it might desirable to early detect
        #        these zero blocks and use a touch[i,j]==false block in ArrayBlock.
        #        How can we detect that we have a zero block?
        block=lazy_map(x->view(x,view_range_i,view_range_j),monolithic_result)
        append!(blocks,[block])
        append!(blocks_coords,[(i,j)])
      end
    end
    cell_grad=lazy_map(BlockMap((nfields,nfields),blocks_coords),blocks...)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

# overloads for AD of SkeletonTriangulation DomainContribution with MultiFields

function FESpaces._change_argument(
  op::typeof(gradient),f,trian::SkeletonTriangulation,uh::MultiFieldFEFunction)

  U = get_fe_space(uh)
  function g(cell_u)
    single_fields_plus = SkeletonCellFieldPair[]
    single_fields_minus = SkeletonCellFieldPair[]
    nfields = length(U.spaces)
    cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
    for i in 1:nfields
      view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
      cell_values_field = lazy_map(a->view(a,view_range),cell_u)
      cf_dual = CellField(U.spaces[i],cell_values_field)
      scfp_plus = SkeletonCellFieldPair(cf_dual, uh[i])
      scfp_minus = SkeletonCellFieldPair(uh[i], cf_dual)
      push!(single_fields_plus,scfp_plus)
      push!(single_fields_minus,scfp_minus)
    end
    xh_plus = MultiFieldCellField(single_fields_plus)
    xh_minus = MultiFieldCellField(single_fields_minus)
    cell_grad_plus = f(xh_plus)
    cell_grad_minus = f(xh_minus)
    get_contribution(cell_grad_plus,trian), get_contribution(cell_grad_minus,trian)
  end
  g
end

function FESpaces._change_argument(
  op::typeof(jacobian),f,trian::SkeletonTriangulation,uh::MultiFieldFEFunction)

  @notimplemented
end

function FESpaces._change_argument(
  op::typeof(hessian),f,trian::SkeletonTriangulation,uh::MultiFieldFEFunction)

  @notimplemented
end

function _restructure_cell_grad!(
  cell_grad, uh::MultiFieldFEFunction, trian::SkeletonTriangulation)

  monolithic_result=cell_grad
  blocks = [] # TO-DO type unstable. How can I infer the type of its entries?
  nfields = length(uh.fe_space.spaces)
  cell_dofs_field_offsets=_get_cell_dofs_field_offsets(uh)
  for i in 1:nfields
    view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
    block_plus = lazy_map(x->view(x[1],view_range),monolithic_result)
    block_minus = lazy_map(x->view(x[2],view_range),monolithic_result)
    block = lazy_map(BlockMap(2,[1,2]),block_plus,block_minus)
    push!(blocks,block)
  end
  cell_grad=lazy_map(BlockMap(nfields,collect(1:nfields)),blocks...)
  cell_grad
end
