
"""
    get_metric_data(trian::Triangulation)

The metric of `trian` as the arguments the metric aware operators take, ready to
be splatted into them:

    hodge_star(ω, get_metric_data(Ω)...)
"""
function get_metric_data(trian::Triangulation)
  model = get_background_model(trian)
  _metric_data(trian, model, get_cell_metric(model))
end

_metric_data(trian::Triangulation, model, ::UniformScaling) = (LinearAlgebra.I,)

function _metric_data(trian::Triangulation, model, model_cell_g::AbstractArray)
  _ensures_full_rank_metric(trian, model)

  mtrian = Triangulation(model)
  model_cell_g_inv = get_cell_inv_metric(model)
  x = testitem(testitem(get_cell_coordinates(mtrian)))
  _check_metric_variance(Cov, return_type(testitem(model_cell_g),     x), "get_cell_metric")
  _check_metric_variance(Con, return_type(testitem(model_cell_g_inv), x), "get_cell_inv_metric")

  m_g     = GenericCellField(model_cell_g,     mtrian, PhysicalDomain())
  m_g_inv = GenericCellField(model_cell_g_inv, mtrian, PhysicalDomain())
  g     = change_domain(m_g    , trian, PhysicalDomain())
  g_inv = change_domain(m_g_inv, trian, PhysicalDomain())

  (g, g_inv, Operation(_sqrt_abs_det)(g))
end

_sqrt_abs_det(g::SymTensorValue) = sqrt(abs(det(g)))

"""
    get_metric_partials_data(trian::Triangulation)

The first order partial derivatives `∂g` of the metric of `trian`, to be fed to

    partials_to_Christoffel_2(get_metric_partials_data(Ω), get_metric_data(Ω)[2])

It is `LinearAlgebra.I` when the metric is the identity, or a cell field.
"""
function get_metric_partials_data(trian::Triangulation)
  model = get_background_model(trian)
  _metric_partials_data(trian, model, get_cell_metric_partials(model))
end

_metric_partials_data(trian::Triangulation, model, ::UniformScaling) = LinearAlgebra.I

function _metric_partials_data(trian::Triangulation, model, model_cell_∂g::AbstractArray)
  _ensures_full_rank_metric(trian, model)
  m_∂g = GenericCellField(model_cell_∂g, Triangulation(model), PhysicalDomain())
  change_domain(m_∂g, trian, PhysicalDomain())
end

# A triangulation of lower dimension than its model, or a model whose cells are
# embedded in a space of higher dimension, carries a metric *induced* on it, which
# is not the cell metric, so both are rejected unless the metric is the identity.
function _ensures_full_rank_metric(trian::Triangulation, model)
  Dm = num_cell_dims(model)
  @notimplementedif num_point_dims(model) != Dm """\n
  A metric of a model whose $(Dm)-dimensional cells are embedded in
  $(num_point_dims(model)) dimensions has rank $Dm in the ambient coordinates its
  fields are evaluated at, so `inv` and `det` of it are not its inverse and its
  volume factor. Carrying one takes a model supplied pseudo inverse and pseudo
  determinant, and a measure built from them, and is not implemented.
  """
  @notimplementedif num_cell_dims(trian) != Dm """\n
  The metric of a $(num_cell_dims(trian))-dimensional triangulation of a
  $(Dm)-dimensional model is the metric it induces on its faces, which is not
  the cell metric restricted to them, and is not implemented.
  """
  nothing
end

# The metric aware operators read the variance of g and g_inv rather than their
# order, so a model carries it; an un-varianced dimension is neither.
function _check_metric_variance(::Type{V},::Type{<:SymTensorValue{D}},accessor) where {V,D}
  D isa V && return nothing
  variance = V === Con ? "contravariant" : "covariant"
  throw(ArgumentError(
    "`$accessor` must return $variance tensors, `SymTensorValue{$V(D)}`, got `SymTensorValue{$D}`"))
end

_check_metric_variance(::Type{V},::Type{T},accessor) where {V,T} = throw(ArgumentError(
  "`$accessor` must return `SymTensorValue` metric tensors, got `$T`"))

