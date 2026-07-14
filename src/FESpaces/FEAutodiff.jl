

function gradient(f::Function,uh::FEFunction;kwargs...)
  fuh = f(uh)
  _gradient(f,uh,fuh;kwargs...)
end

function _gradient(f,uh,fuh::AbstractArray;kwargs...)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.

  Make sure that you are using a Measure instead of a CellQuadrature to perform integration.
  """
end

function _gradient(f,uh,fuh::DomainContribution;tag::GridapADTag=get_ad_level(fuh)+1)
  terms = DomainContribution(;ad_level = tag)
  V = get_fe_space(uh)
  T = eltype(get_vector_type(V))
  for trian in get_domains(fuh)
    g = _change_argument(gradient,f,trian,uh)
    cell_u = get_cell_dof_values(uh)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_gradient(T,g,cell_u,cell_id;tag)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function jacobian(f::Function,uh::FEFunction;kwargs...)
  fuh = f(uh)
  _jacobian(f,uh,fuh;kwargs...)
end

function _jacobian(f,uh,fuh::AbstractArray;kwargs...)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.

  Make sure that you are using a Measure instead of a CellQuadrature to perform integration.
  """
end

function _jacobian(f,uh,fuh::DomainContribution;tag::GridapADTag=get_ad_level(fuh)+1)
  terms = DomainContribution(;ad_level=tag)
  V = get_fe_space(uh)
  T = eltype(get_vector_type(V))
  for trian in get_domains(fuh)
    g = _change_argument(jacobian,f,trian,uh)
    cell_u = get_cell_dof_values(uh)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_jacobian(T,g,cell_u,cell_id;tag)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

"""
    hessian(f::Function, uh::FEFunction)
"""
function hessian(f::Function,uh::FEFunction;kwargs...)
  fuh = f(uh)
  _hessian(f,uh,fuh;kwargs...)
end

function _hessian(f,uh,fuh::AbstractArray;kwargs...)
  @unreachable """\n
  In order to perform AD on a Function taking a FEFunction as argument, such Function
  has to return a DomainContribution.

  Make sure that you are using a Measure instead of a CellQuadrature to perform integration.
  """
end

function _hessian(f,uh,fuh::DomainContribution;tag::GridapADTag=get_ad_level(fuh)+1)
  terms = DomainContribution(;ad_level = tag+1) # Two levels consumed
  V = get_fe_space(uh)
  T = eltype(get_vector_type(V))
  for trian in get_domains(fuh)
    g = _change_argument(hessian,f,trian,uh)
    cell_u = get_cell_dof_values(uh)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_hessian(T,g,cell_u,cell_id;tag)
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

  # Note: In the case where `strian` does not fully cover `ttrian`,
  # tface_to_sface will have negative indices.
  # The negative indices will be dealt with within `autodiff_array_reindex`
  k = 1
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  tface_to_sface = zeros(Int32,length(tface_to_mface))
  for (tface,mface) in enumerate(tface_to_mface)
    sface = mface_to_sface[mface]
    if sface > 0
      tface_to_sface[tface] = sface
    else
      tface_to_sface[tface] = -k
      k += 1
    end
  end

  return tface_to_sface
end

# Skeleton AD

function _change_argument(op,f,trian::SkeletonTriangulation,uh)
  U = get_fe_space(uh)
  function g(cell_u)
    uh_dual = CellField(U,cell_u)
    cf_plus = SkeletonCellFieldPair(uh_dual,uh)
    cf_minus = SkeletonCellFieldPair(uh,uh_dual)
    cg_plus = f(cf_plus)
    cg_minus = f(cf_minus)
    get_contribution(cg_plus,trian), get_contribution(cg_minus,trian)
  end
  g
end

function _compute_cell_ids(uh,ttrian::SkeletonTriangulation)
  plus  = _compute_cell_ids(uh,ttrian.plus)
  minus = _compute_cell_ids(uh,ttrian.minus)
  SkeletonPair(plus,minus)
end

# We collect the derivatives for the plus and minus sides separately,
# which returns ydual_θ = df/duᶿ for θ ∈ {+, -}
# We them merge them into a 2-block BlockVector, so that we obtain
#   result = [df/du⁺, df/du⁻]
function Arrays.autodiff_array_gradient(V, a, i_to_x, j_to_i::SkeletonPair; tag=default_tag(ForwardDiff.gradient,a))
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)

  # dual output of both sides at once
  j_to_ydual_plus, j_to_ydual_minus = a(i_to_xdual)

  # Work for plus side
  j_to_cfg_plus = Arrays.autodiff_array_reindex(i_to_cfg,j_to_i.plus)
  j_to_result_plus = lazy_map(AutoDiffMap(),j_to_cfg_plus,j_to_ydual_plus)

  # Work for minus side
  j_to_cfg_minus = Arrays.autodiff_array_reindex(i_to_cfg,j_to_i.minus)
  j_to_result_minus = lazy_map(AutoDiffMap(),j_to_cfg_minus,j_to_ydual_minus)

  return _skeleton_autodiff_merge_gradient(j_to_result_plus,j_to_result_minus)
end

function _skeleton_autodiff_merge_gradient(j_to_result_plus,j_to_result_minus)
  # Assemble on SkeletonTriangulation expects an array of interior of facets
  # where each entry is a 2-block BlockVector with the first block being the
  # contribution of the plus side and the second, the one of the minus side
  is_single_field = eltype(eltype(j_to_result_plus)) <: Number
  k = is_single_field ? BlockMap(2,[1,2]) : Arrays.BlockBroadcasting(BlockMap(2,[1,2]))
  lazy_map(k,j_to_result_plus,j_to_result_minus)
end

# We collect the derivatives for the plus and minus sides separately,
# which returns ydual_θ = [dr⁺/duᶿ, dr⁻/duᶿ] for θ ∈ {+, -}
# We them merge them as columns into a 2x2 block matrix, so that we obtain
# ydual = [dr⁺/du⁺ dr⁺/du⁻] = [ydual_plus, ydual_minus]
#         [dr⁻/du⁺ dr⁻/du⁻]
function Arrays.autodiff_array_jacobian(V, a, i_to_x, j_to_i::SkeletonPair; tag=default_tag(ForwardDiff.jacobian,a))
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian,tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)

  # dual output of both sides at once
  j_to_ydual_plus, j_to_ydual_minus = a(i_to_xdual)

  # Work for plus side
  j_to_cfg_plus = Arrays.autodiff_array_reindex(i_to_cfg,j_to_i.plus)
  j_to_result_plus = lazy_map(AutoDiffMap(),j_to_cfg_plus,j_to_ydual_plus)

  # Work for minus side
  j_to_cfg_minus = Arrays.autodiff_array_reindex(i_to_cfg,j_to_i.minus)
  j_to_result_minus = lazy_map(AutoDiffMap(),j_to_cfg_minus,j_to_ydual_minus)

  return _skeleton_autodiff_merge_jacobian(j_to_result_plus,j_to_result_minus)
end

function _skeleton_autodiff_merge_jacobian(j_to_result_plus,j_to_result_minus)
  # Merge the columns into a 2x2 block matrix
  I = [
    [(CartesianIndex(1,), CartesianIndex(1, 1)), (CartesianIndex(2,), CartesianIndex(2, 1))], # Plus  -> First column
    [(CartesianIndex(1,), CartesianIndex(1, 2)), (CartesianIndex(2,), CartesianIndex(2, 2))]  # Minus -> Second column
  ]
  is_single_field = eltype(eltype(j_to_result_plus)) <: AbstractArray
  k = is_single_field ? Arrays.MergeBlockMap((2,2),I) : Arrays.BlockBroadcasting(Arrays.MergeBlockMap((2,2),I))
  lazy_map(k,j_to_result_plus,j_to_result_minus)
end

# Skeleton + Complex
# Evaluate r_dual + im*s, then extract the complex gradient/jacobian
function _skeleton_autodiff_array_complex_work(f, a, i_to_x, j_to_i::SkeletonPair, tag)
  s = lazy_map(Broadcasting(imag),i_to_x)
  r = lazy_map(Broadcasting(real),i_to_x)
  i_to_cfg   = lazy_map(ConfigMap(f,tag),r)
  i_to_rdual = lazy_map(DualizeMap(),i_to_cfg,r)
  j_to_ycdual_plus, j_to_ycdual_minus = a(lazy_map((r,s) -> r + im*s, i_to_rdual, s))
  j_to_cfg_plus  = Arrays.autodiff_array_reindex(i_to_cfg,j_to_i.plus)
  j_to_cfg_minus = Arrays.autodiff_array_reindex(i_to_cfg,j_to_i.minus)
  j_to_result_plus  = lazy_map((u,v) -> u + im*v,
    lazy_map(AutoDiffMap(),j_to_cfg_plus,  lazy_map(Broadcasting(real),j_to_ycdual_plus)),
    lazy_map(AutoDiffMap(),j_to_cfg_plus,  lazy_map(Broadcasting(imag),j_to_ycdual_plus)))
  j_to_result_minus = lazy_map((u,v) -> u + im*v,
    lazy_map(AutoDiffMap(),j_to_cfg_minus, lazy_map(Broadcasting(real),j_to_ycdual_minus)),
    lazy_map(AutoDiffMap(),j_to_cfg_minus, lazy_map(Broadcasting(imag),j_to_ycdual_minus)))

  return j_to_result_plus,j_to_result_minus
end

function Arrays.autodiff_array_gradient(::Type{<:Complex},a,i_to_x,j_to_i::SkeletonPair;tag=default_tag(ForwardDiff.gradient,a))
  j_to_result_plus,j_to_result_minus = _skeleton_autodiff_array_complex_work(ForwardDiff.gradient,a,i_to_x,j_to_i,tag)
  # Merge contributions
  is_single_field = eltype(eltype(j_to_result_plus)) <: Number
  k = is_single_field ? BlockMap(2,[1,2]) : Arrays.BlockBroadcasting(BlockMap(2,[1,2]))
  lazy_map(k,j_to_result_plus,j_to_result_minus)
end

function Arrays.autodiff_array_jacobian(::Type{<:Complex},a,i_to_x,j_to_i::SkeletonPair;tag=default_tag(ForwardDiff.jacobian,a))
  j_to_result_plus,j_to_result_minus = _skeleton_autodiff_array_complex_work(ForwardDiff.jacobian,a,i_to_x,j_to_i,tag)
  # Merge contributions
  I = [
    [(CartesianIndex(1,), CartesianIndex(1, 1)), (CartesianIndex(2,), CartesianIndex(2, 1))], # Plus  -> First column
    [(CartesianIndex(1,), CartesianIndex(1, 2)), (CartesianIndex(2,), CartesianIndex(2, 2))]  # Minus -> Second column
  ]
  is_single_field = eltype(eltype(j_to_result_plus)) <: AbstractArray
  k = is_single_field ? Arrays.MergeBlockMap((2,2),I) : Arrays.BlockBroadcasting(Arrays.MergeBlockMap((2,2),I))
  lazy_map(k,j_to_result_plus,j_to_result_minus)
end