using Gridap


abstract type FEObservationOperator end #<: FEINNType end
abstract type DataMisfitCalculator end #<: FEINNType end

FEObservationOperator(a) = a 

function filter_observation_values(::FEObservationOperator, ::AbstractVector)
  @abstractmethod
end

function _update_free_dof_values!(vh::FEFunction, v::AbstractVector)
  get_free_dof_values(vh) .= v
end

# function _update_free_dof_values!(
#   vh::Union{DistributedCellField,DistributedMultiFieldFEFunction},
#   v::PVector)
#   u = get_free_dof_values(vh)
#   u .= v
#   if u.index_partition !== v.index_partition
#     fetch_vector_ghost_values!(partition(u), map(reverse, u.cache)) |> wait
#   end
# end

"""
An operator that given the DoF values of a FE function returns
a global vector resulting from evaluating the function at a set
of arbitrary observation points within the computational domain.

# Fields
- obs_points: A vector of observation points.
- fe_space: An FE space where points belong to.
- cell_to_points: The mapping from triangulation cells to the points within it.
- cell_w_points: The list of triangulation cells where there is at least one point within it.
- cell_to_weights: FE space shape functions evaluated at observation points, stored cell-wise.
- filtered_indices: The new indices for the points after filtering out invalid points.
- cache: Preallocated vectors for computing the observation values and the adjoints.
"""
struct SequentialFEObservationOperator{T} <: FEObservationOperator
  obs_points::AbstractVector{<:Point}
  fe_space::SingleFieldFESpace
  cell_to_points::Gridap.Arrays.Table
  cells_w_points
  cell_to_weights::AbstractVector{Matrix{T}}
  filtered_indices::Vector{Int32}
  cache
end

function FEObservationOperator(
  points::AbstractVector{<:Point},
  fe_space::SingleFieldFESpace,
  model=nothing)  # model is not necessary
  test_space = fe_space
  if isa(test_space, TrialFESpace)
    test_space = fe_space.space
  end

  Ω = get_triangulation(test_space)
  fidx, point_to_cell = Int32[], Int32[]
  cache = _point_to_cell_cache(KDTreeSearch(), Ω)
  for (i, point) in enumerate(points)
    try
      cid = _point_to_cell!(cache, point)
      push!(fidx, i)
      push!(point_to_cell, cid)
    catch _
      # printstyled("WARNING: "; bold=true, color=:yellow)
      # println("Point $point is not inside any active cell.")
    end
  end
  cell_to_points, _ = make_inverse_table(point_to_cell, num_cells(Ω))
  cells_w_points = findall(i -> (cell_to_points.ptrs[i+1] - cell_to_points.ptrs[i]) != 0,
    1:length(cell_to_points.ptrs)-1)
  cell_to_xs = lazy_map(Broadcasting(Reindex(points[fidx])), cell_to_points)
  cell_map = get_cell_map(Ω)
  inv_map = lazy_map(inverse_map, cell_map)
  cell_xs_at_ref_space = lazy_map(evaluate, inv_map, cell_to_xs)
  cell_basis = get_data(get_fe_basis(test_space))
  cell_to_weights = lazy_map(evaluate, cell_basis, cell_xs_at_ref_space)
  WT = eltype(eltype(cell_to_weights))
  cache = _init_obs_op_cache(WT, fe_space, length(fidx))
  SequentialFEObservationOperator{WT}(
    points, fe_space, cell_to_points,
    cells_w_points, cell_to_weights,
    fidx, cache)
end

function filter_observation_values(
  op::SequentialFEObservationOperator,
  obs_values::AbstractVector)
  @assert length(obs_values) == length(op.obs_points) """
  Number of observation points does NOT match number of observation values!
  """

  obs_values[op.filtered_indices]
end

"""
Given the free DoF values, evaluate the corresponding
FE function at the observation points.

# Arguments
- free_vals: A vector of free DoF values.
"""
function (obsop::SequentialFEObservationOperator)(free_vals)
  obs_vals, uh, _ = obsop.cache
  _update_free_dof_values!(uh, free_vals)
  _evaluate_obs_values!(
    obs_vals,
    obsop.cell_to_points,
    obsop.cells_w_points,
    obsop.cell_to_weights,
    get_cell_dof_values(uh))
  obs_vals
end

function ChainRulesCore.rrule(obsop::SequentialFEObservationOperator, free_vals)
  obs_vals = obsop(free_vals)

  function sequential_obs_operator_pullback(adjoint_parent)
    free_v̄als = obsop.cache[end-1]
    _compute_obs_op_adjoints!(
      obsop.cell_to_weights,
      get_cell_dof_ids(obsop.fe_space),
      obsop.cell_to_points,
      obsop.cells_w_points,
      adjoint_parent,
      free_v̄als)
    return NoTangent(), free_v̄als
  end

  return obs_vals, sequential_obs_operator_pullback
end

"""
Given the free and Dirichlet DoF values, evaluate the
corresponding FE function at the observation points.

# Arguments
- free_vals: A vector of the DoF values.
- diri_vals: A vector of the Dirichlet values.
"""
function (op::SequentialFEObservationOperator)(free_vals, diri_vals)
  obs_vals, uh, _ = op.cache
  @show free_vals
  @show length(uh.free_values)
  _update_free_dof_values!(uh, free_vals)
  uh.dirichlet_values .= diri_vals
  _evaluate_obs_values!(
    obs_vals,
    op.cell_to_points,
    op.cells_w_points,
    op.cell_to_weights,
    get_cell_dof_values(uh))
  obs_vals
end

function ChainRulesCore.rrule(obsop::SequentialFEObservationOperator, free_vals, diri_vals)
  obs_vals = obsop(free_vals, diri_vals)

  function sequential_obs_operator_pullback(adjoint_parent)
    free_v̄als, diri_v̄als = obsop.cache[end-1], obsop.cache[end]
    _compute_obs_op_adjoints!(
      obsop.cell_to_weights,
      get_cell_dof_ids(obsop.fe_space),
      obsop.cell_to_points,
      obsop.cells_w_points,
      adjoint_parent,
      free_v̄als,
      diri_v̄als)
    return NoTangent(), free_v̄als, diri_v̄als
  end

  return obs_vals, sequential_obs_operator_pullback
end

# struct DistributedFEObservationOperator <: FEObservationOperator
#   obs_points::AbstractVector{<:Point}
#   fe_space::DistributedFESpace
#   model::DistributedDiscreteModel
#   operators::AbstractVector{SequentialFEObservationOperator}
#   cache
# end

# function FEObservationOperator(
#   points::AbstractVector{<:Point},
#   fe_space::DistributedFESpace,
#   model::DistributedDiscreteModel)

#   obs_ops = map(space -> FEObservationOperator(points, space), fe_space.spaces)
#   pobs_idx = _compute_obs_data_indices(model, obs_ops, length(points))
#   pobs_vals = PVector(map(op -> op.cache[1], obs_ops), pobs_idx, nothing)
#   pfree_v̄als = PVector(map(op -> op.cache[end-1], obs_ops), partition(fe_space.gids))
#   diri_spaces = map(s -> DirichletFESpace(s), fe_space.spaces)
#   diri_gids = generate_gids(model, diri_spaces)
#   pdiri_v̄als = nothing
#   if length(diri_gids) > 1
#     pdiri_v̄als = PVector(map(op -> op.cache[end], obs_ops), partition(diri_gids))
#   end
#   cache = (pobs_vals, pfree_v̄als, pdiri_v̄als)
#   DistributedFEObservationOperator(points, fe_space, model, obs_ops, cache)
# end

# function filter_observation_values(op::DistributedFEObservationOperator, obs_values::AbstractVector)
#   @assert length(obs_values) == length(op.obs_points) """
#   Number of observation points does NOT match number of observation values!
#   """

#   filtered_values = map(op -> obs_values[op.filtered_indices], op.operators)
#   PVector(filtered_values, op.cache[1].index_partition, nothing)
# end

# function (obsop::DistributedFEObservationOperator)(free_vals)
#   obs_vals = obsop.cache[1]
#   map(obsop.operators, partition(free_vals)) do op, fv
#     op(fv)
#   end
#   obs_vals
# end

# function ChainRulesCore.rrule(obsop::DistributedFEObservationOperator, free_vals)
#   obs_vals = obsop(free_vals)

#   function distributed_obs_operator_pullback(adjoint_parents)
#     pfree_v̄als = obsop.cache[end-1]
#     map(obsop.operators, partition(adjoint_parents)) do op, adj_parent
#       free_v̄als = op.cache[end-1]
#       _compute_obs_op_adjoints!(
#         op.cell_to_weights,
#         get_cell_dof_ids(op.fe_space),
#         op.cell_to_points,
#         op.cells_w_points,
#         adj_parent,
#         free_v̄als)
#     end
#     return NoTangent(), pfree_v̄als
#   end

#   return obs_vals, distributed_obs_operator_pullback
# end

# function (obsop::DistributedFEObservationOperator)(free_vals, diri_vals)
#   obs_vals = obsop.cache[1]
#   map(obsop.operators, partition(free_vals), partition(diri_vals)) do op, fv, dv
#     op(fv, dv)
#   end
#   obs_vals
# end

# function ChainRulesCore.rrule(obsop::DistributedFEObservationOperator, free_vals, diri_vals)
#   obs_vals = obsop(free_vals, diri_vals)

#   function distributed_obs_operator_pullback(adj_parents)
#     pfree_v̄als, pdiri_v̄als = obsop.cache[end-1], obsop.cache[end]
#     map(obsop.operators, partition(adj_parents)) do op, adj_parent
#       free_v̄als, diri_v̄als = op.cache[end-1], op.cache[end]
#       _compute_obs_op_adjoints!(
#         op.cell_to_weights,
#         get_cell_dof_ids(op.fe_space),
#         op.cell_to_points,
#         op.cells_w_points,
#         adj_parent,
#         free_v̄als,
#         diri_v̄als)
#     end
#     return NoTangent(), pfree_v̄als, pdiri_v̄als
#   end

#   return obs_vals, distributed_obs_operator_pullback
# end


# struct MultiFieldDataMisfitCalculator <: DataMisfitCalculator
#   fe_space::Union{DistributedMultiFieldFESpace,MultiFieldFESpace}
#   obs_operators::Vector{<:FEObservationOperator}
#   obs_values::Vector{<:AbstractVector}
#   loss_functions::Vector{<:Function}
#   weights::Vector{<:Real}
#   cache
# end

# function MultiFieldDataMisfitCalculator(
#   fe_space::Union{DistributedMultiFieldFESpace,MultiFieldFESpace},
#   obs_operators::Vector{<:FEObservationOperator},
#   obs_values::Vector{<:AbstractVector};
#   loss_functions::Vector{<:Function}=Function[], weights=Float64[])
#   nops = length(obs_operators)
#   isempty(weights) && (weights = ones(nops))
#   isempty(loss_functions) && (loss_functions = [norm2 for _ in 1:nops])
#   _validate_multi_field_data_misfit_cal_inputs(fe_space, obs_operators, obs_values, loss_functions, weights)

#   supported_losses = Set([norm1, norm2, norm_sqr])
#   @assert all([lf in supported_losses for lf in loss_functions]) """Supported
#   loss functions are `norm1`, `norm2`, and `norm_sqr` from `LinearAlgebra`."""

#   cache = _init_multi_field_data_misfit_cal_cache(fe_space, obs_values)
#   MultiFieldDataMisfitCalculator(
#     fe_space, obs_operators, obs_values,
#     loss_functions, weights, cache)
# end

# function (cal::MultiFieldDataMisfitCalculator)(u::AbstractVector{T}) where {T<:Real}
#   upartitions, ustart_idx, _ = cal.cache
#   total_loss = zero(T)
#   for (i, loss_func) in enumerate(cal.loss_functions)
#     wi, opi, ui, ovi = cal.weights[i], cal.obs_operators[i], upartitions[i], cal.obs_values[i]
#     _update_subvector_from_vector!(ui, u, ustart_idx, i; transform=partition)
#     total_loss += wi * loss_func(opi(ui) - ovi)
#   end
#   total_loss
# end

# function ChainRulesCore.rrule(cal::MultiFieldDataMisfitCalculator, u::PVector)
#   upartitions, ustart_idx, ūall, errs = cal.cache
#   total_loss, num_field = zero(eltype(u)), length(cal.obs_operators)
#   opsbacks, errbacks = [], Vector{Function}(undef, num_field)
#   for (i, upart) in enumerate(upartitions)
#     _update_subvector_from_vector!(upart, u, ustart_idx, i; transform=partition)
#     ops, ovi = cal.obs_operators[i].operators, cal.obs_values[i]
#     opsback = map(ops, partition(ovi), partition(upart), partition(errs[i])) do op, ov, up, e
#       opv, opback = ChainRulesCore.rrule(op, up)
#       e .= opv - ov
#       opback
#     end
#     push!(opsbacks, opsback)
#     ji, errback = ChainRulesCore.rrule(cal.loss_functions[i], errs[i])
#     total_loss += cal.weights[i] * ji
#     errbacks[i] = errback
#   end

#   function data_misfit_calculator_pullback(j̄)
#     wone = one(eltype(cal.weights))
#     for (i, opsback) in enumerate(opsbacks)
#       ērr = errbacks[i](j̄)[2]
#       map(opsback, partition(ūall), ustart_idx, partition(ērr)) do opback, ūa, uidx, ē
#         cal.weights[i] != wone && (ē .*= cal.weights[i])
#         ūpart = opback(ē)[2]
#         _update_vector_from_subvector!(ūa, ūpart, uidx, i)
#       end
#     end
#     NoTangent(), ūall
#   end

#   total_loss, data_misfit_calculator_pullback
# end

# function ChainRulesCore.rrule(cal::MultiFieldDataMisfitCalculator, u::Vector{T}) where {T}
#   upartitions, ustart_idx, ūall, errs = cal.cache
#   nparts, total_loss = length(upartitions), zero(T)
#   errbacks, opbacks = Vector{Function}(undef, nparts), Vector{Function}(undef, nparts)
#   for (i, upart) in enumerate(upartitions)
#     _update_subvector_from_vector!(upart, u, ustart_idx, i)
#     opv, opback = ChainRulesCore.rrule(cal.obs_operators[i], upart)
#     opbacks[i] = opback
#     errs[i] .= opv - cal.obs_values[i]
#     ji, errback = ChainRulesCore.rrule(cal.loss_functions[i], errs[i])
#     total_loss += cal.weights[i] * ji
#     errbacks[i] = errback
#   end

#   function data_misfit_calculator_pullback(j̄)
#     wone = one(eltype(cal.weights))
#     for (i, opback) in enumerate(opbacks)
#       ērr = unthunk(errbacks[i](j̄)[2])
#       cal.weights[i] != wone && (ērr .*= cal.weights[i])
#       ūpart = opback(ērr)[2]
#       _update_vector_from_subvector!(ūall, ūpart, ustart_idx, i)
#     end
#     NoTangent(), ūall
#   end

#   total_loss, data_misfit_calculator_pullback
# end


"""
Init the cache for the observation operator.

# Arguments
- WT: the type of the observation values.
- fe_space: the FE space.
- nobs: number of observations.

# returns
- obs_vals: a preallocated vector to store observation values.
- uh: an FE function defined on the FE space.
- free_v̄als: a preallocated vector to store the adjoint of the free dof values.
- diri_v̄als: a preallocated vector to store the adjoint of the Dirichlet dof values.
"""
function _init_obs_op_cache(WT, fe_space, nobs)
  obs_vals = zeros(WT, nobs)
  VT = eltype(get_vector_type(fe_space))
  free_v̄als = Vector{VT}(undef, num_free_dofs(fe_space))
  diri_v̄als = Vector{VT}(undef, num_dirichlet_dofs(fe_space))
  uh = FEFunction(fe_space, similar(free_v̄als), copy(get_dirichlet_dof_values(fe_space)))
  (obs_vals, uh, free_v̄als, diri_v̄als)
end

function _evaluate_obs_values!(
  obs_vals,
  cell_to_points::Gridap.Arrays.Table,
  cells_w_points,
  cell_to_weights,
  cell_to_dof_values)
  cell_evals = lazy_map(*, cell_to_weights, cell_to_dof_values)
  Gridap.FESpaces._free_and_dirichlet_values_fill!(
    obs_vals,
    nothing, # No Dirichlet values in the present scenario
    array_cache(cell_evals),
    array_cache(cell_to_points),
    cell_evals,
    cell_to_points,
    cells_w_points)
end

function _free_and_dirichlet_values_accum!(
  free_vals,
  cache_vals,
  cache_dofs,
  cell_vals,
  cell_dofs,
  cells,
  diri_vals=nothing)

  free_vals .= zero(eltype(free_vals))
  compute_diri_vals = diri_vals !== nothing
  if compute_diri_vals
    diri_vals .= zero(eltype(diri_vals))
  end

  for cell in cells
    vals = getindex!(cache_vals, cell_vals, cell)
    dofs = getindex!(cache_dofs, cell_dofs, cell)
    for (i, dof) in enumerate(dofs)
      val = vals[i]
      if dof > 0
        free_vals[dof] += val
      elseif compute_diri_vals && dof < 0
        diri_vals[-dof] += val
      elseif dof == 0
        @unreachable "dof ids either positive or negative, not zero"
      end
    end
  end

end

function _compute_obs_op_adjoints!(
  cell_to_weights,
  cell_to_dof_ids,
  cell_to_points,
  cells_w_points,
  adjoint_parent,
  free_v̄als,
  diri_v̄als=nothing)

  cell_vals = lazy_map(transpose, lazy_map(Broadcasting(Reindex(adjoint_parent)), cell_to_points))
  cell_vals_mul_weights = lazy_map(*, cell_vals, cell_to_weights)
  _free_and_dirichlet_values_accum!(
    free_v̄als,
    array_cache(cell_vals_mul_weights),
    array_cache(cell_to_dof_ids),
    cell_vals_mul_weights,
    cell_to_dof_ids,
    cells_w_points,
    diri_v̄als)
end

function _compute_obs_data_indices(model, obs_ops, npoint)
  cell_gids = partition(get_cell_gids(model))
  map(cell_gids, obs_ops) do cgid, op
    l2o = local_to_owner(cgid)
    owners = Vector{Int32}(undef, length(op.filtered_indices))
    for (i, pts) in enumerate(op.cell_to_points)
      for ipt in pts
        owners[ipt] = l2o[i]
      end
    end
    LocalIndices(npoint, first(own_to_owner(cgid)), Int.(op.filtered_indices), owners)
  end
end

# function _validate_multi_field_data_misfit_cal_inputs(
#   fe_space::DistributedMultiFieldFESpace,
#   obs_operators,
#   obs_values,
#   loss_functions,
#   weights)
#   @assert length(fe_space.field_fe_space) ==
#           length(obs_operators) ==
#           length(obs_values) ==
#           length(loss_functions) ==
#           length(weights) """Number of single-field FE spaces, number of
#           observation operators and values, number of weights,
#           number of loss functions must be the same!"""
#   @assert all([fe_space.field_fe_space[i].gids === obs_operators[i].fe_space.gids
#                for i in eachindex(obs_operators)]) """The FE spaces in the
#                multi-field FE space must match those in the observation operators"""
# end

function _validate_multi_field_data_misfit_cal_inputs(
  fe_space::MultiFieldFESpace,
  obs_operators,
  obs_values,
  loss_functions,
  weights)
  @assert length(fe_space.spaces) ==
          length(obs_operators) ==
          length(obs_values) ==
          length(loss_functions) ==
          length(weights) """Number of single-field FE spaces, number of
          observation operators and values, number of weights,
          number of loss functions must be the same!"""
  @assert all([fe_space.spaces[i] === obs_operators[i].fe_space
               for i in eachindex(obs_operators)]) """The FE spaces in the
               multi-field FE space must match those in the observation operators"""
end

function _init_multi_field_data_misfit_cal_cache(fe_space::MultiFieldFESpace, obs_values)
  PT, nspaces = eltype(get_vector_type(fe_space)), length(fe_space.spaces)
  ustart_idx, upartitions = Vector{Int32}(undef, nspaces), Vector{Vector{PT}}(undef, nspaces)
  uidx = 1
  for i in 1:nspaces
    ustart_idx[i] = uidx
    nu = num_free_dofs(fe_space.spaces[i])
    upartitions[i] = Vector{PT}(undef, nu)
    uidx += nu
  end
  ūall = vcat(upartitions...)
  errs = map(ov -> similar(ov), obs_values)
  (upartitions, ustart_idx, ūall, errs)
end

# function _init_multi_field_data_misfit_cal_cache(fe_space::DistributedMultiFieldFESpace, obs_values)
#   PT = eltype(get_vector_type(fe_space))
#   nspaces = length(fe_space.field_fe_space)
#   ustart_idx = map(_ -> [Int32(1)], fe_space.part_fe_space)
#   upartitions = Vector{PVector{Vector{PT}}}(undef, nspaces)
#   for (i, space) in enumerate(fe_space.field_fe_space)
#     upar_vec = map(partition(space.gids), ustart_idx) do gids, uidx
#       nlocal = length(gids)
#       i < nspaces && push!(uidx, uidx[end] + nlocal)
#       Vector{PT}(undef, nlocal)
#     end
#     upartitions[i] = PVector(upar_vec, partition(space.gids))
#   end

#   ūall_vec = map(ps -> Vector{PT}(undef, num_free_dofs(ps)), fe_space.part_fe_space)
#   ūall = PVector(ūall_vec, partition(fe_space.gids))
#   errs = map(ov -> similar(ov), obs_values)
#   (upartitions, ustart_idx, ūall, errs)
# end