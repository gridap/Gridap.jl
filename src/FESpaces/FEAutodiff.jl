

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
    cell_u = get_cell_dof_values(uh)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_gradient(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

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
    cell_u = get_cell_dof_values(uh)
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
    cell_u = get_cell_dof_values(uh)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_hessian(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function _change_argument(op,f,trian,uh::SingleFieldFEFunction)
  U = get_fe_space(uh)
  function g(cell_u)
    cf = CellField(U,cell_u)
    cell_grad = f(cf)
    get_contribution(cell_grad,trian)
  end
  g
end

#= AD for DomainContribution involving SkeletonTriangulation

- Following are the constructs for performing gradient of DomainContribution
  involving integrations over SkeletonTriangulation (Λ)
- The current approach followed to achieve the above is performing the Gridap
  way of AD for plus and minus sides of the FEFunction occuring at Λ separately,
  and combining the result. So as to Dualize only either plus side or minus
  side of CellField/FEFunction we introduce the SkeletonCellFieldPair, which
  stores two CellFields, one of which in the use case here is the dualized
  version of the other.
- Currently, Jacobian and hence Hessian are not yet fully implemented and work-
  in-progress.
- AD for integration over SkeletonTriangulation has not yet been implemented
  for the MultiField case

=#

function _change_argument(
  op,f,
  trian::SkeletonTriangulation,
  uh::SingleFieldFEFunction)

  U = get_fe_space(uh)
  function g(cell_u)
    uh_dual = CellField(U,cell_u)
    scfp_plus = SkeletonCellFieldPair(uh_dual, uh, trian)
    scfp_minus = SkeletonCellFieldPair(uh, uh_dual, trian)
    cell_grad_plus = f(scfp_plus)
    cell_grad_minus = f(scfp_minus)
    get_contribution(cell_grad_plus,trian), get_contribution(cell_grad_minus,trian)
  end
  g
end

# In the earlier version, SkeletonPair was returned. Although this is a perfect
# fit for this situation, the dispatch based on SkeletonPair is not possible as
# it is not imported into src/Arrays/ due to circula dependency situation
function _compute_cell_ids(uh,ttrian::SkeletonTriangulation)
  tcells_plus  = _compute_cell_ids(uh,ttrian.plus)
  tcells_minus = _compute_cell_ids(uh,ttrian.minus)
  # SkeletonPair(tcells_plus,tcells_minus)
  (tcells_plus, tcells_minus)
end
