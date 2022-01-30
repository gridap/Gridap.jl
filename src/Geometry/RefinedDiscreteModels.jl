# This implementation of the newest vertex bisection algorithm is based on the
# following article:

# LONG CHEN (2008). SHORT IMPLEMENTATION OF BISECTION IN MATLAB.
# In Recent Advances in Computational Sciences. WORLD SCIENTIFIC.
# DOI: 10.1142/9789812792389_0020
using Gridap.Arrays
using SparseArrays
using TimerOutputs
using Random

# Create a TimerOutput, this is the main type that keeps track of everything.
const to = TimerOutput()

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

function shift_to_first(v::AbstractArray{Ti}, i::Ti) where {Ti<:Integer}
  circshift(v, -(i - 1))
end

function sort_longest_edge!(
    elem::Table{Ti},
    node::AbstractArray{<:VectorValue},
    NT::Ti,
) where {Ti<:Integer}
    edgelength = zeros(NT, 3)
    #node_v = [[v[1], v[2]] for v in node]
    #node_v = vcat(node_v'...)
    for t = 1:NT
      elem_t = elem[t][:]
        for (j, e) in enumerate(elem_t)
            arr = filter(x -> x != e, elem_t)
            diff = sqrt(sum((node[arr[1]] - node[arr[2]]) .^ 2))
            #diff = sqrt(sum((node_v[arr[1], :] - node_v[arr[2], :]).^2))
            edgelength[t, j] = diff
        end
    end
    max_indices = findmax(edgelength, dims = 2)[2]
    for t = 1:NT
      shifted = shift_to_first(elem[t][:], Ti(max_indices[t][2]))
      for j = 1:3
        elem.data[elem.ptrs[t] + j - 1] = shifted[j]
      end
    end
end

function setup_markers_and_nodes!(
    node::AbstractArray{<:VectorValue},
    elem::Table{Ti},
    d2p::SparseMatrixCSC{Ti,Ti},
    dualedge::SparseMatrixCSC{Ti},
    NE::Ti,
    η_arr::AbstractArray{<:AbstractFloat},
    θ::AbstractFloat,
  ) where {Ti <: Integer}
  total_η = sum(η_arr)
  partial_η = 0
  sorted_η_idxs = sortperm(-η_arr)
  marker = zeros(Ti, NE)
  # Loop over global triangle indices
  for t = sorted_η_idxs
    if (partial_η > θ * total_η)
      break
    end
    need_to_mark = true
    # Get triangle index with next largest error
    #ct = sorted_η_idxs[t]
    while (need_to_mark)
      # Base point
      base = d2p[elem[t][2], elem[t][3]]
      # Already marked
      if marker[base] > 0
        need_to_mark = false
      else
        # Get the estimator contribution for the current triangle
        partial_η = partial_η + η_arr[t]
        # Increase the number of nodes to add new midpoint
        N = size(node, 1) + 1
        # The marker of the current elements is this node
        marker[d2p[elem[t][2], elem[t][3]]] = N
        # Coordinates of new node
        elem2 = elem[t][2]
        elem3 = elem[t][3]
        midpoint = get_midpoint(node[[elem2, elem3], :])
        node = push!(node, midpoint)
        t = dualedge[elem[t][3], elem[t][2]]
        # There is no dual edge here, go to next triangle index
        if t == 0
          need_to_mark = false
        end
      end
    end
  end
  node, marker
end

function divide!(elem::AbstractArray, t::Ti, p::AbstractArray{Ti}) where {Ti <: Integer}
  new_row = [p[4], p[3], p[1]]
  update_row = [p[4], p[1], p[2]]
  push!(elem, new_row)
  elem[t] = update_row
  #elem = append_tables_globally(elem, Table([new_row]))
  #for j = 1:3
  #  elem.data[elem.ptrs[t] - 1 + j] = update_row[j]
  #end
  #elem[t][:] = [p[4] p[1] p[2]]
  elem
end

function bisect(
    d2p::SparseMatrixCSC{Ti,Ti},
    elem::AbstractArray,
    marker::AbstractArray{Ti},
    NT::Ti,
  ) where {Ti<:Integer}
  for t = UnitRange{Ti}(1:NT)
    base = d2p[elem[t][2], elem[t][3]]
    if (marker[base] > 0)
      p = vcat(elem[t][:], marker[base])
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

function build_edges(elem::Table{Ti}) where {Ti <: Integer}
  edge = zeros(Ti, 3*length(elem), 2)
  for t = 1:length(elem)
    off = 3*(t-1)
    edge[off+1,:] = [elem[t][1] elem[t][2]]
    edge[off+2, :] = [elem[t][1] elem[t][3]]
    edge[off+3, :] = [elem[t][2] elem[t][3]]
  end
  unique(sort!(edge, dims = 2), dims = 1)
end

function build_directed_dualedge(elem::Table{Ti}, N::Ti, NT::Ti) where {Ti <: Integer}
  dualedge = spzeros(Ti, N, N)
  for t = 1:NT
    dualedge[elem[t][1], elem[t][2]] = t
    dualedge[elem[t][2], elem[t][3]] = t
    dualedge[elem[t][3], elem[t][1]] = t
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

function test_against_top(face::Table{<:Integer}, top::GridTopology, d::Integer)
  face_vec = sort.([face[i][:] for i = 1:size(face, 1)])
  face_top = sort.(get_faces(top, d, 0))
  issetequal_bitvec = issetequal(face_vec, face_top)
  @assert all(issetequal_bitvec)
end

function test_against_top(face::Matrix{<:Integer}, top::GridTopology, d::Integer)
  face_vec = sort.([face[i, :] for i = 1:size(face, 1)])
  face_top = sort.(get_faces(top, d, 0))
  issetequal_bitvec = issetequal(face_vec, face_top)
  #@assert all(issetequal_bitvec)
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

function sort_ccw(cell_coords::AbstractArray{<:VectorValue})
  midpoint = get_midpoint(cell_coords)
  offset_coords = cell_coords .- midpoint
  sortperm(offset_coords, by = atan)
end

function sort_cell_node_ids_ccw(
    cell_node_ids::Table{<:Integer},
    node_coords::AbstractArray{<:VectorValue},
  )
  #cell_node_ids_ccw = vcat(cell_node_ids'...)
  cell_node_ids_ccw = cell_node_ids
  #@show cell_node_ids_ccw
  for (i, cell) in enumerate(cell_node_ids)
    cell_coords = node_coords[cell]
    perm = sort_ccw(cell_coords)
    #cell_node_ids_ccw[i][:] = cell[perm]
    permed = cell[perm]
    for j = 1:3
      cell_node_ids_ccw.data[cell_node_ids_ccw.ptrs[i] + j - 1] = permed[j]
    end
  end
  cell_node_ids_ccw
end

"""
node_coords == node, cell_node_ids == elem in Long Chen's notation
"""
function newest_vertex_bisection(
  top::GridTopology,
  node_coords::Vector,
  cell_node_ids::Table{Ti},
  η_arr::AbstractArray{<:AbstractFloat},
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
  #test_against_top(elem, top, 2)
  @timeit to "edges" edge = build_edges(elem)
  NE::Ti = size(edge, 1)
  @timeit to "dualedge" dualedge = build_directed_dualedge(elem, N, NT)
  d2p = dual_to_primal(edge, NE, N)
  # Make sure edge is consistent with GridTopology
  #@timeit to "test" test_against_top(edge, top, 1)
  @timeit to "markers" node_coords, marker =
    setup_markers_and_nodes!(node_coords, elem, d2p, dualedge, NE, η_arr, θ)
  @timeit to "copy elem" elem = Vector{Vector}(elem)
  @timeit to "bisect" cell_node_ids = bisect(d2p, elem, marker, NT)
  #@show cell_node_ids
  @timeit to "recreate" cell_node_ids = Table([c for c in cell_node_ids])
  node_coords, cell_node_ids
end

# step 1
function newest_vertex_bisection(
  grid::Grid,
  top::GridTopology,
  η_arr::AbstractArray{<:AbstractFloat},
  θ::AbstractFloat,
  sort_flag::Bool,
)
  node_coords = get_node_coordinates(grid)
  cell_node_ids = get_cell_node_ids(grid)
  if sort_flag
    cell_node_ids_ccw = sort_cell_node_ids_ccw(cell_node_ids, node_coords)
  else
    cell_node_ids_ccw = cell_node_ids
  end
  @timeit to "newest" node_coords_ref, cell_node_ids_ref =
    newest_vertex_bisection(top, node_coords, cell_node_ids_ccw, η_arr, θ, sort_flag)
  # TODO: Should not convert to matrix and back to Table
  #cell_node_ids_ref = Table([c for c in eachrow(cell_node_ids_ref)])
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
  η_arr::AbstractArray{<:AbstractFloat};
  θ = 1.0, # corresponds to uniform refinement
  sort_flag = false,
)
  # Not sure if necessary to keep old model unchanged. For my tests I use this
  model_c = deepcopy(model)
  grid = get_grid(model_c)
  top = get_grid_topology(model_c)
  ref_grid = newest_vertex_bisection(grid, top, η_arr, θ, sort_flag)
  ref_topo = GridTopology(ref_grid)
  #ref_labels = # Compute them from the original labels (This is perhaps the most tedious part)
  ref_labels = FaceLabeling(ref_topo)
  DiscreteModel(ref_grid, ref_topo, ref_labels)
end

function build_refined_models(
  model::DiscreteModel,
  Nsteps::Integer,
  θ::AbstractFloat,
  est::Estimator,
)
  reset_timer!(to)
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
  #@show to
  model_refs
end
