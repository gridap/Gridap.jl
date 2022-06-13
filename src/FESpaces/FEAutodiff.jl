

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

function _compute_cell_ids(uh,ttrian::SkeletonTriangulation)
  tcells_plus  = _compute_cell_ids(uh,ttrian.plus)
  tcells_minus = _compute_cell_ids(uh,ttrian.minus)
  SkeletonPair(tcells_plus,tcells_minus)
end

## overloads for AD of SkeletonTriangulation DomainContribution ##

#= Notes regarding the placement of below AD functions for Skeleton Integration

- Earlier, the autodiff_array_### family of functions have been placed in the
  `src/Arrays/Autodiff.jl`, but as below we are leveraging `BlockMap` for the
  construction derivatives of Skeleton integration terms, this cannot be
  placed in `Arrays/Autodiff.jl` due to circular dependency, this is due to the
  fact that `BlockMap` belongs to Gridap.Fields which uses/imports functions
  and constructs from Gridap.Arrays, so using `BlockMap` in Gridap.Arrays would
  create a circular dependency.
- So as to have everything working we will need `use` some of the constructs
  from Gridap.Arrays, `Gridap.Geometry.SkeletonPair` and
  `Gridap.CellData.SkeletonCellFieldPair`
- This comes with an added advantage of being able use `SkeletonPair` in
  `_compute_cell_ids` for `SkeletonTriangulation`
- we need to also think if all of AD can be moved into a separate module in
  Gridap, where we can use or import required functionalities from other
  modules without any circular dependency
=#
function autodiff_array_gradient(
  a, i_to_x, j_to_i::SkeletonPair)

  i_to_xdual = lazy_map(DualizeMap(ForwardDiff.gradient),i_to_x)

  # dual output of both sides at once
  j_to_ydual_plus, j_to_ydual_minus = a(i_to_xdual)

  # Work for plus side
  j_to_x_plus = lazy_map(Reindex(i_to_x),j_to_i.plus)
  j_to_cfg_plus = lazy_map(ConfigMap(ForwardDiff.gradient),j_to_x_plus)
  j_to_result_plus = lazy_map(AutoDiffMap(ForwardDiff.gradient),
                              j_to_ydual_plus,j_to_x_plus,j_to_cfg_plus)

  # Work for minus side
  j_to_x_minus = lazy_map(Reindex(i_to_x),j_to_i.minus)
  j_to_cfg_minus = lazy_map(ConfigMap(ForwardDiff.gradient),j_to_x_minus)
  j_to_result_minus = lazy_map(AutoDiffMap(ForwardDiff.gradient),
                               j_to_ydual_minus,j_to_x_minus,j_to_cfg_minus)

  # Assemble on SkeletonTriangulation expects an array of interior of facets
  # where each entry is a 2-block BlockVector with the first block being the
  # contribution of the plus side and the second, the one of the minus side
  lazy_map(BlockMap(2,[1,2]),j_to_result_plus,j_to_result_minus)
end

function autodiff_array_jacobian(
  a,i_to_x,j_to_i::SkeletonPair)
  @notimplemented
end

function autodiff_array_hessian(
  a,i_to_x,i_to_j::SkeletonPair)
  @notimplemented
end
