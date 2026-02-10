function _update_free_dof_values!(
  vh::Union{DistributedCellField,DistributedMultiFieldFEFunction},
  v::PVector)
  u = get_free_dof_values(vh)
  u .= v
  if u.index_partition !== v.index_partition
    fetch_vector_ghost_values!(partition(u), map(reverse, u.cache)) |> wait
  end
end

struct DistributedFEObservationOperator <: FEObservationOperator
  obs_points::AbstractVector{<:Point}
  fe_space::DistributedFESpace
  model::DistributedDiscreteModel
  operators::AbstractVector{SequentialFEObservationOperator}
  cache
end

function FEObservationOperator(
  points::AbstractVector{<:Point},
  fe_space::DistributedFESpace,
  model::DistributedDiscreteModel)

  obs_ops = map(space -> FEObservationOperator(points, space), fe_space.spaces)
  pobs_idx = _compute_obs_data_indices(model, obs_ops, length(points))
  pobs_vals = PVector(map(op -> op.cache[1], obs_ops), pobs_idx, nothing)
  pfree_v̄als = PVector(map(op -> op.cache[end-1], obs_ops), partition(fe_space.gids))
  diri_spaces = map(s -> DirichletFESpace(s), fe_space.spaces)
  diri_gids = generate_gids(model, diri_spaces)
  pdiri_v̄als = nothing
  if length(diri_gids) > 1
    pdiri_v̄als = PVector(map(op -> op.cache[end], obs_ops), partition(diri_gids))
  end
  cache = (pobs_vals, pfree_v̄als, pdiri_v̄als)
  DistributedFEObservationOperator(points, fe_space, model, obs_ops, cache)
end

function filter_observation_values(op::DistributedFEObservationOperator, obs_values::AbstractVector)
  @assert length(obs_values) == length(op.obs_points) """
  Number of observation points does NOT match number of observation values!
  """

  filtered_values = map(op -> obs_values[op.filtered_indices], op.operators)
  PVector(filtered_values, op.cache[1].index_partition, nothing)
end

function (obsop::DistributedFEObservationOperator)(free_vals)
  obs_vals = obsop.cache[1]
  map(obsop.operators, partition(free_vals)) do op, fv
    op(fv)
  end
  obs_vals
end

function ChainRulesCore.rrule(obsop::DistributedFEObservationOperator, free_vals)
  obs_vals = obsop(free_vals)

  function distributed_obs_operator_pullback(adjoint_parents)
    pfree_v̄als = obsop.cache[end-1]
    map(obsop.operators, partition(adjoint_parents)) do op, adj_parent
      free_v̄als = op.cache[end-1]
      _compute_obs_op_adjoints!(
        op.cell_to_weights,
        get_cell_dof_ids(op.fe_space),
        op.cell_to_points,
        op.cells_w_points,
        adj_parent,
        free_v̄als)
    end
    return NoTangent(), pfree_v̄als
  end

  return obs_vals, distributed_obs_operator_pullback
end

function (obsop::DistributedFEObservationOperator)(free_vals, diri_vals)
  obs_vals = obsop.cache[1]
  map(obsop.operators, partition(free_vals), partition(diri_vals)) do op, fv, dv
    op(fv, dv)
  end
  obs_vals
end

function ChainRulesCore.rrule(obsop::DistributedFEObservationOperator, free_vals, diri_vals)
  obs_vals = obsop(free_vals, diri_vals)

  function distributed_obs_operator_pullback(adj_parents)
    pfree_v̄als, pdiri_v̄als = obsop.cache[end-1], obsop.cache[end]
    map(obsop.operators, partition(adj_parents)) do op, adj_parent
      free_v̄als, diri_v̄als = op.cache[end-1], op.cache[end]
      _compute_obs_op_adjoints!(
        op.cell_to_weights,
        get_cell_dof_ids(op.fe_space),
        op.cell_to_points,
        op.cells_w_points,
        adj_parent,
        free_v̄als,
        diri_v̄als)
    end
    return NoTangent(), pfree_v̄als, pdiri_v̄als
  end

  return obs_vals, distributed_obs_operator_pullback
end


struct MultiFieldDataMisfitCalculator <: DataMisfitCalculator
  fe_space::Union{DistributedMultiFieldFESpace,MultiFieldFESpace}
  obs_operators::Vector{<:FEObservationOperator}
  obs_values::Vector{<:AbstractVector}
  loss_functions::Vector{<:Function}
  weights::Vector{<:Real}
  cache
end

function MultiFieldDataMisfitCalculator(
  fe_space::Union{DistributedMultiFieldFESpace,MultiFieldFESpace},
  obs_operators::Vector{<:FEObservationOperator},
  obs_values::Vector{<:AbstractVector};
  loss_functions::Vector{<:Function}=Function[], weights=Float64[])
  nops = length(obs_operators)
  isempty(weights) && (weights = ones(nops))
  isempty(loss_functions) && (loss_functions = [norm2 for _ in 1:nops])
  _validate_multi_field_data_misfit_cal_inputs(fe_space, obs_operators, obs_values, loss_functions, weights)

  supported_losses = Set([norm1, norm2, norm_sqr])
  @assert all([lf in supported_losses for lf in loss_functions]) """Supported
  loss functions are `norm1`, `norm2`, and `norm_sqr` from `LinearAlgebra`."""

  cache = _init_multi_field_data_misfit_cal_cache(fe_space, obs_values)
  MultiFieldDataMisfitCalculator(
    fe_space, obs_operators, obs_values,
    loss_functions, weights, cache)
end

function (cal::MultiFieldDataMisfitCalculator)(u::AbstractVector{T}) where {T<:Real}
  upartitions, ustart_idx, _ = cal.cache
  total_loss = zero(T)
  for (i, loss_func) in enumerate(cal.loss_functions)
    wi, opi, ui, ovi = cal.weights[i], cal.obs_operators[i], upartitions[i], cal.obs_values[i]
    _update_subvector_from_vector!(ui, u, ustart_idx, i; transform=partition)
    total_loss += wi * loss_func(opi(ui) - ovi)
  end
  total_loss
end

function ChainRulesCore.rrule(cal::MultiFieldDataMisfitCalculator, u::PVector)
  upartitions, ustart_idx, ūall, errs = cal.cache
  total_loss, num_field = zero(eltype(u)), length(cal.obs_operators)
  opsbacks, errbacks = [], Vector{Function}(undef, num_field)
  for (i, upart) in enumerate(upartitions)
    _update_subvector_from_vector!(upart, u, ustart_idx, i; transform=partition)
    ops, ovi = cal.obs_operators[i].operators, cal.obs_values[i]
    opsback = map(ops, partition(ovi), partition(upart), partition(errs[i])) do op, ov, up, e
      opv, opback = ChainRulesCore.rrule(op, up)
      e .= opv - ov
      opback
    end
    push!(opsbacks, opsback)
    ji, errback = ChainRulesCore.rrule(cal.loss_functions[i], errs[i])
    total_loss += cal.weights[i] * ji
    errbacks[i] = errback
  end

  function data_misfit_calculator_pullback(j̄)
    wone = one(eltype(cal.weights))
    for (i, opsback) in enumerate(opsbacks)
      ērr = errbacks[i](j̄)[2]
      map(opsback, partition(ūall), ustart_idx, partition(ērr)) do opback, ūa, uidx, ē
        cal.weights[i] != wone && (ē .*= cal.weights[i])
        ūpart = opback(ē)[2]
        _update_vector_from_subvector!(ūa, ūpart, uidx, i)
      end
    end
    NoTangent(), ūall
  end

  total_loss, data_misfit_calculator_pullback
end

function ChainRulesCore.rrule(cal::MultiFieldDataMisfitCalculator, u::Vector{T}) where {T}
  upartitions, ustart_idx, ūall, errs = cal.cache
  nparts, total_loss = length(upartitions), zero(T)
  errbacks, opbacks = Vector{Function}(undef, nparts), Vector{Function}(undef, nparts)
  for (i, upart) in enumerate(upartitions)
    _update_subvector_from_vector!(upart, u, ustart_idx, i)
    opv, opback = ChainRulesCore.rrule(cal.obs_operators[i], upart)
    opbacks[i] = opback
    errs[i] .= opv - cal.obs_values[i]
    ji, errback = ChainRulesCore.rrule(cal.loss_functions[i], errs[i])
    total_loss += cal.weights[i] * ji
    errbacks[i] = errback
  end

  function data_misfit_calculator_pullback(j̄)
    wone = one(eltype(cal.weights))
    for (i, opback) in enumerate(opbacks)
      ērr = unthunk(errbacks[i](j̄)[2])
      cal.weights[i] != wone && (ērr .*= cal.weights[i])
      ūpart = opback(ērr)[2]
      _update_vector_from_subvector!(ūall, ūpart, ustart_idx, i)
    end
    NoTangent(), ūall
  end

  total_loss, data_misfit_calculator_pullback
end

function _validate_multi_field_data_misfit_cal_inputs(
  fe_space::DistributedMultiFieldFESpace,
  obs_operators,
  obs_values,
  loss_functions,
  weights)
  @assert length(fe_space.field_fe_space) ==
          length(obs_operators) ==
          length(obs_values) ==
          length(loss_functions) ==
          length(weights) """Number of single-field FE spaces, number of
          observation operators and values, number of weights,
          number of loss functions must be the same!"""
  @assert all([fe_space.field_fe_space[i].gids === obs_operators[i].fe_space.gids
               for i in eachindex(obs_operators)]) """The FE spaces in the
               multi-field FE space must match those in the observation operators"""
end

function _init_multi_field_data_misfit_cal_cache(fe_space::DistributedMultiFieldFESpace, obs_values)
  PT = eltype(get_vector_type(fe_space))
  nspaces = length(fe_space.field_fe_space)
  ustart_idx = map(_ -> [Int32(1)], fe_space.part_fe_space)
  upartitions = Vector{PVector{Vector{PT}}}(undef, nspaces)
  for (i, space) in enumerate(fe_space.field_fe_space)
    upar_vec = map(partition(space.gids), ustart_idx) do gids, uidx
      nlocal = length(gids)
      i < nspaces && push!(uidx, uidx[end] + nlocal)
      Vector{PT}(undef, nlocal)
    end
    upartitions[i] = PVector(upar_vec, partition(space.gids))
  end

  ūall_vec = map(ps -> Vector{PT}(undef, num_free_dofs(ps)), fe_space.part_fe_space)
  ūall = PVector(ūall_vec, partition(fe_space.gids))
  errs = map(ov -> similar(ov), obs_values)
  (upartitions, ustart_idx, ūall, errs)
end