

function gradient(f::Function,uh::FEFunction)
  fuh = f(uh)
  _gradient(f,uh,fuh)
end

function _gradient(f,uh,fuh::AbstractArray)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.

  Make sure that you are using a Measure instead of a CellQuadrature to perform integration.
  """
end

function _gradient(f,uh,fuh::DomainContribution)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = _change_argument(gradient,f,trian,uh)
    cell_u = _get_cell_dof_values(uh,trian)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_gradient(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function jacobian(f::Function,uh::FEFunction)
  fuh = f(uh)
  _jacobian(f,uh,fuh)
end

function _jacobian(f,uh,fuh::AbstractArray)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.

  Make sure that you are using a Measure instead of a CellQuadrature to perform integration.
  """
end

function _jacobian(f,uh,fuh::DomainContribution)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = _change_argument(jacobian,f,trian,uh)
    cell_u = _get_cell_dof_values(uh,trian)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_jacobian(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function hessian(f::Function,uh::FEFunction)
  fuh = f(uh)
  _hessian(f,uh,fuh)
end

function _hessian(f,uh,fuh::AbstractArray)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.

  Make sure that you are using a Measure instead of a CellQuadrature to perform integration.
  """
end

function _hessian(f,uh,fuh::DomainContribution)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = _change_argument(hessian,f,trian,uh)
    cell_u = _get_cell_dof_values(uh,trian)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_hessian(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function _change_argument(op,f,trian,uh)
  U = get_fe_space(uh)
  function g(cell_u)
    cf = CellField(U,cell_u)
    cell_grad = f(cf)
    get_contribution(cell_grad,trian)
  end
  g
end

_get_cell_dof_values(uh,trian) = get_cell_dof_values(uh)

function _compute_cell_ids(uh,ttrian)
  strian = get_triangulation(uh)
  if strian === ttrian
    return collect(IdentityVector(Int32(num_cells(strian))))
  end
  @check is_change_possible(strian,ttrian)
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  @notimplementedif !isa(sglue,FaceToFaceGlue)
  @notimplementedif !isa(tglue,FaceToFaceGlue)
  scells = IdentityVector(Int32(num_cells(strian)))
  mcells = extend(scells,sglue.mface_to_tface)
  tcells = lazy_map(Reindex(mcells),tglue.tface_to_mface)
  collect(tcells)
end

# Skeleton AD

function _change_argument(op,f,trian::SkeletonTriangulation,uh)
  U = get_fe_space(uh)
  function g(cell_u)
    cf = SkeletonCellFieldPair(U,cell_u)
    cell_grad = f(cf)
    get_contribution(cell_grad,trian)
  end
  g
end

function _get_cell_dof_values(uh,trian::SkeletonTriangulation)
  plus  = _get_cell_dof_values(uh,trian.plus)
  minus = _get_cell_dof_values(uh,trian.minus)
  lazy_map(BlockMap(2,[1,2]),plus,minus)
end

function _compute_cell_ids(uh,ttrian::SkeletonTriangulation)
  plus  = _compute_cell_ids(uh,ttrian.plus)
  minus = _compute_cell_ids(uh,ttrian.minus)
  lazy_map(BlockMap(2,[1,2]),plus,minus)
end

function CellData.SkeletonCellFieldPair(V::FESpace,cell_values)
  plus = CellField(V,lazy_map(x->x.array[1],cell_values))
  minus = CellField(V,lazy_map(x->x.array[2],cell_values))
  CellData.SkeletonCellFieldPair(plus,minus)
end

function Arrays.lazy_map(r::Reindex, ids::LazyArray{<:Fill{BlockMap{1}}})
  k = ids.maps.value
  plus = lazy_map(r,ids.args[1])
  # minus = lazy_map(r,ids.args[2])
  # lazy_map(k,plus,minus)
  plus
end
