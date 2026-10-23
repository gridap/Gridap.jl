
"""
    get_metric_data(trian::Triangulation)

The metric of `trian` as the arguments the metric aware operators take, ready to
be splatted into them:

    hodge_star(ω, get_metric_data(Ω)...)

Either `(LinearAlgebra.I,)`, the identity metric, or `(g, g_inv, _sqrt_abs_det)`.

A triangulation whose cells do not fill the ambient space carries the metric its
background *induces* on them and not the cell metric of that background. Which
of the two applies is read off the dimensions, see [`_is_embedded`](@ref): every
triangulation of the tree reads its cells in the frames of its background model,
so none of them needs a method of its own. What they do not share is the
orientation their cells carry, [`_induced_orientation`](@ref), which has no
default and is where a new triangulation declares itself.
"""
function get_metric_data(trian::Triangulation)
  model = get_background_model(trian)
  _is_embedded(trian) &&
    return _induced_metric_data(trian, model, get_cell_metric(model))
  _metric_data(trian, model, get_cell_metric(model))
end

# A triangulation of lower dimension than its background carries the metric the
# inclusion pulls back from it, not that background's cell metric: that is what
# its fields are paired and starred with. A surface mesh, whose own cells are
# embedded in a higher dimensional space, is deliberately not this case -- its
# fields have ambient components, which the ambient identity pairs correctly and
# whose area element the integration map supplies. What frame a triangulation's
# components live in is not a function of its dimensions alone.
_is_embedded(trian::Triangulation, model) =
  num_point_dims(model) == num_cell_dims(model) &&
  num_cell_dims(trian)  <  num_cell_dims(model)

_is_embedded(trian::Triangulation) = _is_embedded(trian, get_background_model(trian))

# `g_Γ = ι*g = Jι ⋅ g ⋅ Jιᵀ`, in the reference cell coordinates of `trian`,
# `Jι` being the differential of the inclusion read in those coordinates and
# the ambient frame, see [`_get_inclusion_jac`](@ref).
#
# Slot 3 is the component of the volume form of `g_Γ` in those coordinates, not
# merely its volume factor: `g_Γ` is a Gram matrix, so `det g_Γ > 0` and the
# magnitude is derivable from `g_Γ`, while the sign is not. A `⋆` takes an
# orientation as well as a metric, and the one a face carries is the one its
# background induces on it, so slot 3 is signed by `ε`. That is what makes the
# two writings of a top form integral agree,
#
#   ∫( ω )dΓ = ∫( ⋆ω )dΓ,
#
# the first carrying `ε` in the weights, the second in the `⋆`. Its magnitude is
# the area element `meas(Jι)` that the integration map already applies, so it is
# not added to a scalar integrand a second time, see [`_metric_weights`](@ref).
function _induced_metric_data(trian::Triangulation, model, ::UniformScaling)
  g = Operation(_induced_metric)(_get_inclusion_jac(trian))
  ε = CellField(_induced_orientation(trian), trian)
  (g, Operation(inv)(g), ε * Operation(_sqrt_abs_det)(g))
end

function _induced_metric_data(trian::Triangulation, model, cell_g)
  @notimplemented """\n
  The metric a $(num_cell_dims(trian))-dimensional triangulation inherits from
  its $(num_cell_dims(model))-dimensional background is the pullback of the
  background metric by the inclusion. For a non identity background metric its
  volume factor is not the `meas` of the cell map that the integration map
  applies, and is not implemented.
  """
end

# `ι*g` is the congruence of `g` by the differential of the inclusion, the
# pullback of a twice covariant tensor, which for the Kronecker delta of an
# ambient frame is [`congruent_prod`](@ref) of that differential alone.
_induced_metric(Jι::MultiValue, g::MultiValue) = pullback(g, Jι)
_induced_metric(Jι::MultiValue) = congruent_prod(LinearAlgebra.I, Jι)

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

"""
    get_metric_partials_data(trian::Triangulation)

The first order partial derivatives `∂g` of the metric of `trian`, to be fed to

    partials_to_Christoffel_2(get_metric_partials_data(Ω), get_metric_data(Ω)[2])

It is `LinearAlgebra.I` when the metric is the identity, or a cell field.

They are the derivatives of the metric [`get_metric_data`](@ref) returns, and
like it they are a property of the background model rather than of the
triangulation, except that an embedded triangulation has none: the derivatives
of an induced metric are not those of the background's.
"""
function get_metric_partials_data(trian::Triangulation)
  model = get_background_model(trian)
  @notimplementedif _is_embedded(trian) """\n
  The metric of a $(num_cell_dims(trian))-dimensional triangulation embedded in
  $(num_point_dims(trian)) dimensions is the one its background induces on it,
  see [`get_metric_data`](@ref), so its partial derivatives are those of that
  induced metric and not those of the background's, which are not implemented.
  The operators that need them, `δ` in particular, are not available there.
  """
  _metric_partials_data(trian, model, get_cell_metric_partials(model))
end

_metric_partials_data(trian::Triangulation, model, ::UniformScaling) = LinearAlgebra.I

function _metric_partials_data(trian::Triangulation, model, model_cell_∂g::AbstractArray)
  _ensures_full_rank_metric(trian, model)
  m_∂g = GenericCellField(model_cell_∂g, Triangulation(model), PhysicalDomain())
  change_domain(m_∂g, trian, PhysicalDomain())
end

# A triangulation of lower dimension than its model, or a model whose cells are
# embedded in a higher dimensional space, carries a metric *induced* on it and
# not the cell metric, so both are rejected unless that metric is the identity.
# [`_is_embedded`](@ref) catches them upstream: this is the safety net for a
# `_metric_data` reached by another path.
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

