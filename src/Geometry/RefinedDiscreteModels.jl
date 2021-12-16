using Gridap.Arrays
using SparseArrays
using Random

# For testing only
abstract type Estimator end

struct ConstantEst <: Estimator
  val::Float64
end

struct RandomEst <: Estimator
  function RandomEst(seed)
    Random.seed!(seed)
    new()
  end
end

# For testing only
compute_estimator(est::RandomEst, ncells) = rand(ncells)
compute_estimator(est::ConstantEst, ncells) = fill(est.val, ncells)

function shift_to_first(v::Vector{Ti}, i::Ti) where {Ti<:Integer}
  circshift(v, -(i - 1))
end

function sort_longest_edge!(
    elem::Matrix{Ti},
    node::Vector{<:VectorValue},
    NT::Ti,
) where {Ti<:Integer}
    edgelength = zeros(NT, 3)
    #node_v = [[v[1], v[2]] for v in node]
    #node_v = vcat(node_v'...)
    for i = 1:NT
        elem_i = elem[i, :]
        for (j, e) in enumerate(elem_i)
            arr = filter(x -> x != e, elem_i)
            diff = sqrt(sum((node[arr[1]] - node[arr[2]]) .^ 2))
            #diff = sqrt(sum((node_v[arr[1], :] - node_v[arr[2], :]).^2))
            edgelength[i, j] = diff
        end
    end
    max_indices = findmax(edgelength, dims = 2)[2]
    for i = 1:NT
      elem[i, :] = shift_to_first(elem[i, :], Ti(max_indices[i][2]))
    end
end

function setup_markers_and_nodes!(
    node::Vector{<:VectorValue},
    elem::Matrix{Ti},
    d2p::SparseMatrixCSC{Ti,Ti},
    dualedge::SparseMatrixCSC{Ti},
    NE::Ti,
    η_arr::Vector{<:AbstractFloat},
    θ::AbstractFloat,
  ) where {Ti <: Integer}
  total_η = sum(η_arr)
  partial_η = 0
  sorted_η_idxs = sortperm(-η_arr)
  marker = zeros(Ti, NE)
  # Loop over global triangle indices
  for t_idx = sorted_η_idxs
    if (partial_η > θ * total_η)
      break
    end
    need_to_mark = true
    # Get triangle index with next largest error
    #ct = sorted_η_idxs[t]
    while (need_to_mark)
      # Base point
      base = d2p[elem[t_idx, 2], elem[t_idx, 3]]
      # Already marked
      if marker[base] > 0
        need_to_mark = false
      else
        # Get the estimator contribution for the current triangle
        partial_η = partial_η + η_arr[t_idx]
        # Increase the number of nodes to add new midpoint
        N = size(node, 1) + 1
        # The marker of the current elements is this node
        marker[d2p[elem[t_idx, 2], elem[t_idx, 3]]] = N
        # Coordinates of new node
        midpoint = get_midpoint(node[elem[t_idx, [2 3], :]])
        node = push!(node, midpoint)
        t_idx = dualedge[elem[t_idx, 3], elem[t_idx, 2]]
        # There is no dual edge here, go to next triangle index
        if t_idx == 0
          need_to_mark = false
        end
      end
    end
  end
  node, marker
end

function divide!(elem::Matrix{Ti}, t::Ti, p::Vector{Ti}) where {Ti <: Integer}
  elem = [elem; [p[4] p[3] p[1]]]
  elem[t, :] = [p[4] p[1] p[2]]
  elem
end

function bisect(
    d2p::SparseMatrixCSC{Ti,Ti},
    elem::Matrix{Ti},
    marker::Vector{Ti},
    NT::Ti,
  ) where {Ti<:Integer}
  for t = UnitRange{Ti}(1:NT)
    base = d2p[elem[t, 2], elem[t, 3]]
    if (marker[base] > 0)
      p = vcat(elem[t, :], marker[base])
      elem = divide!(elem, t, p)
      left = d2p[p[1], p[2]]
      right = d2p[p[3], p[1]]
      if (marker[right] > 0)
        cur_size::Ti = size(elem, 1)
        elem = divide!(elem, cur_size, [p[4], p[3], p[1], marker[right]])
      end
      if (marker[left] > 0)
        elem = divide!(elem, t, [p[4], p[1], p[2], marker[left]])
      end
    end
  end
  elem
end

function build_edges(elem::Matrix{<:Integer})
  edge = [elem[:, [1, 2]]; elem[:, [1, 3]]; elem[:, [2, 3]]]
  unique(sort!(edge, dims = 2), dims = 1)
end

function build_directed_dualedge(elem::Matrix{Ti}, N::Ti, NT::Ti) where {Ti <: Integer}
  dualedge = spzeros(Ti, N, N)
  for t = 1:NT
    dualedge[elem[t, 1], elem[t, 2]] = t
    dualedge[elem[t, 2], elem[t, 3]] = t
    dualedge[elem[t, 3], elem[t, 1]] = t
  end
  dualedge
end

function dual_to_primal(edge::Matrix{Ti}, NE::Ti, N::Ti) where {Ti <: Integer}
  d2p = spzeros(Ti, Ti, N, N)
  for k = 1:NE
    i = edge[k, 1]
    j = edge[k, 2]
    d2p[i, j] = k
    d2p[j, i] = k
  end
  d2p
end

function test_against_top(face::Matrix{<:Integer}, top::GridTopology, d::Integer)
  face_vec = sort.([face[i, :] for i = 1:size(face, 1)])
  face_top = sort.(get_faces(top, d, 0))
  issetequal_bitvec = issetequal(face_vec, face_top)
  @assert all(issetequal_bitvec)
end

function get_midpoint(ngon::AbstractArray{<:VectorValue})
  sum(ngon) / length(ngon)
end

function Base.atan(v::VectorValue{2,T}) where {T<:AbstractFloat}
  atan(v[2], v[1])
end

function Base.:^(v::VectorValue{N,T}, r::Integer) where {T, N}
  VectorValue([v[i]^r for i = 1:N]...)
end

function sort_ccw(cell_coords::Vector{<:VectorValue})
  midpoint = get_midpoint(cell_coords)
  offset_coords = cell_coords .- midpoint
  sortperm(offset_coords, by = atan)
end

function sort_cell_node_ids_ccw(
    cell_node_ids::Table{<:Integer},
    node_coords::Vector{<:VectorValue},
  )
  cell_node_ids_ccw = vcat(cell_node_ids'...)
  #cell_node_ids_ccw = cell_node_ids
  #@show cell_node_ids_ccw
  for (i, cell) in enumerate(cell_node_ids)
    cell_coords = node_coords[cell]
    perm = sort_ccw(cell_coords)
    cell_node_ids_ccw[i, :] = cell[perm]
  end
  cell_node_ids_ccw
end

"""
node_coords == node, cell_node_ids == elem in Long Chen's notation
"""
function newest_vertex_bisection(
  top::GridTopology,
  node_coords::Vector,
  cell_node_ids::Matrix{Ti},
  η_arr::Vector{<:AbstractFloat},
  θ::AbstractFloat,
  sort_flag::Bool,
) where {Ti <: Integer}
  # Number of nodes
  N::Ti = size(node_coords, 1)
  elem = cell_node_ids
  # Number of cells (triangles in 2D)
  NT::Ti = size(elem, 1)
  @assert length(η_arr) == NT
  if sort_flag
    sort_longest_edge!(elem, node_coords, NT)
  end
  # Make sure elem is consistent with GridTopology
  test_against_top(elem, top, 2)
  edge = build_edges(elem)
  NE::Ti = size(edge, 1)
  dualedge = build_directed_dualedge(elem, N, NT)
  d2p = dual_to_primal(edge, NE, N)
  # Make sure edge is consistent with GridTopology
  test_against_top(edge, top, 1)
  node_coords, marker =
    setup_markers_and_nodes!(node_coords, elem, d2p, dualedge, NE, η_arr, θ)
  cell_node_ids = bisect(d2p, elem, marker, NT)
  node_coords, cell_node_ids
end

# step 1
function newest_vertex_bisection(
  grid::Grid,
  top::GridTopology,
  η_arr::AbstractVector{<:AbstractFloat},
  θ::AbstractFloat,
  sort_flag::Bool,
)
  node_coords = get_node_coordinates(grid)
  cell_node_ids = get_cell_node_ids(grid)
  if sort_flag
    cell_node_ids_ccw = sort_cell_node_ids_ccw(cell_node_ids, node_coords)
  else
    cell_node_ids_ccw = vcat(cell_node_ids'...)
  end
  node_coords_ref, cell_node_ids_ref =
    newest_vertex_bisection(top, node_coords, cell_node_ids_ccw, η_arr, θ, sort_flag)
  # TODO: Should not convert to matrix and back to Table
  cell_node_ids_ref = Table([c for c in eachrow(cell_node_ids_ref)])
  reffes = get_reffes(grid)
  cell_types = get_cell_type(grid)
  # TODO : Gracefully handle cell_types
  new_cell_types = fill(1, length(cell_node_ids_ref) - length(cell_node_ids))
  append!(cell_types, new_cell_types)
  UnstructuredGrid(node_coords_ref, cell_node_ids_ref, reffes, cell_types)
end

# step 2
function newest_vertex_bisection(
  model::DiscreteModel,
  η_arr::AbstractVector{<:AbstractFloat};
  θ = 1.0, # corresponds to uniform refinement
  sort_flag = false,
)
  grid = get_grid(model)
  top = get_grid_topology(model)
  ref_grid = newest_vertex_bisection(grid, top, η_arr, θ, sort_flag)
  ref_topo = GridTopology(ref_grid)
  #ref_labels = # Compute them from the original labels (This is perhaps the most tedious part)
  ref_labels = FaceLabeling(ref_topo)
  DiscreteModel(ref_grid, ref_topo, ref_labels)
end

function build_refined_models(
  domain::Tuple,
  partition::Tuple,
  Nsteps::Integer,
  θ::AbstractFloat,
  est::Estimator,
)
  model = CartesianDiscreteModel(domain, partition)
  model = simplexify(model)
  model_refs = Vector{DiscreteModel}(undef, Nsteps)
  cell_map = get_cell_map(get_triangulation(model))
  ncells = length(cell_map)
  η_arr = compute_estimator(est, ncells)
  model_refs[1] = newest_vertex_bisection(model, η_arr; sort_flag = true, θ = θ)
  for i = 1:(Nsteps - 1)
    cell_map = get_cell_map(get_triangulation(model_refs[i]))
    ncells = length(cell_map)
    η_arr = compute_estimator(est, ncells)
    model_refs[i + 1] =
      newest_vertex_bisection(model_refs[i], η_arr; sort_flag = false, θ = θ)
  end
  model_refs
end
