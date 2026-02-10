abstract type FEObservationOperator end 

function filter_observation_values(::FEObservationOperator, ::AbstractVector)
  @abstractmethod
end

function _update_free_dof_values!(vh::FEFunction, v::AbstractVector)
  get_free_dof_values(vh) .= v
end

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
  cell_to_points::Table
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
  cell_to_points::Table,
  cells_w_points,
  cell_to_weights,
  cell_to_dof_values)
  cell_evals = lazy_map(*, cell_to_weights, cell_to_dof_values)
  _free_and_dirichlet_values_fill!(
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