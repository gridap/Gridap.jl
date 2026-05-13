# ===================================================================
# Abstract interface
# ===================================================================

"""
    reindex_free_and_dirichlet_dof_ids(space::FESpace, free_dof_ids, dir_dof_ids)

Return a new `FESpace` obtained by reindexing the free and Dirichlet DOF IDs.
`free_dof_ids[old_free] = new_free` and `dir_dof_ids` uses the same convention
(negative values preserve the Dirichlet sign encoding of `cell_dof_ids`).
Concrete implementations exist for `UnconstrainedFESpace` and `PolytopalFESpace`.
"""
function reindex_free_and_dirichlet_dof_ids(space::FESpace,free_dof_ids,dir_dof_ids)
  @abstractmethod
end

"""
    reindex_free_dof_ids(space::FESpace, free_dof_ids)

Apply a given permutation of the free DOFs, leaving Dirichlet DOF IDs unchanged.
`free_dof_ids[old_free] = new_free`.
"""
function reindex_free_dof_ids(space::FESpace, free_dof_ids)
  reindex_free_and_dirichlet_dof_ids(space, free_dof_ids, -1:-1:-num_dirichlet_dofs(space))
end

# TODO: deprecate on next major release
function renumber_free_and_dirichlet_dof_ids(space::FESpace,args...)
  reindex_free_and_dirichlet_dof_ids(space,args...)
end

# ===================================================================
# George-Liu pseudo-peripheral node finder (shared by RCM and Sloan).
#
# Starts from `seed` and iterates BFS from the minimum-degree node in
# the last level set until the eccentricity stops increasing.
# Returns the pseudo-peripheral node (the previous BFS root).
# ===================================================================
function _george_liu_peripheral!(bfs!::F, queue, level, degrees, seed) where {F}
  n, last = bfs!(seed)
  ecc     = level[queue[n]]
  current = seed
  while true
    best, best_deg = 0, typemax(Int)
    for qi in last:n
      v = queue[qi]
      degrees[v] < best_deg && (best_deg = degrees[v]; best = v)
    end
    n2, last2 = bfs!(best)
    ecc2 = level[queue[n2]]
    ecc2 <= ecc && return current
    ecc = ecc2; n = n2; last = last2; current = best
  end
end

# ===================================================================
# Reverse Cuthill-McKee permutation.
#
# Uses George-Liu for pseudo-peripheral node selection. Sorts each
# adjacency row by ascending degree on a private copy so CM visits
# lower-degree neighbours first without mutating the caller's adj.
#
# State encoding: 0 = unvisited, cg = visited in generation cg, -1 = done.
# Hot-loop condition state[u] % UInt < cg % UInt covers all three cases.
# ===================================================================
function _rcm_perm(adj::Table)
  nfree   = length(adj)
  degrees = diff(adj.ptrs)

  sorted_data = copy(adj.data)
  adj = Table(sorted_data, adj.ptrs)
  for i in 1:nfree
    sort!(view(sorted_data, datarange(adj, i)), by=j -> degrees[j])
  end

  queue = Vector{Int}(undef, nfree)
  level = Vector{Int}(undef, nfree)
  state = zeros(Int, nfree)
  gen   = Ref(0)

  function bfs!(start)
    cg = (gen[] += 1)
    state[start] = cg; level[start] = 0
    queue[1] = start; qhead = qtail = 1
    max_level = 0
    while qhead <= qtail
      v = queue[qhead]; qhead += 1; lv = level[v]
      for k in datarange(adj, v)
        u = adj.data[k]
        if state[u] % UInt < cg % UInt
          state[u] = cg; level[u] = lv + 1
          lv + 1 > max_level && (max_level = lv + 1)
          qtail += 1; queue[qtail] = u
        end
      end
    end
    last = qtail
    while last > 1 && level[queue[last-1]] == max_level; last -= 1; end
    return qtail, last
  end

  perm = Vector{Int}(undef, nfree)
  idx  = 1

  while idx <= nfree
    seed, seed_deg = 0, typemax(Int)
    for i in 1:nfree
      state[i] == 0 && degrees[i] < seed_deg && (seed_deg = degrees[i]; seed = i)
    end
    peripheral = _george_liu_peripheral!(bfs!, queue, level, degrees, seed)
    n, _ = bfs!(peripheral)
    for qi in 1:n
      v = queue[qi]; state[v] = -1; perm[idx] = v; idx += 1
    end
  end

  reverse!(perm)
  iperm = Vector{Int}(undef, nfree)
  for i in 1:nfree; iperm[perm[i]] = i; end
  return iperm
end

# ===================================================================
# Sloan's algorithm.
#
# Minimises the active front (wavefront) rather than bandwidth.
# Uses a priority queue with f(i) = eff_w1·dist_E(i) + eff_w2·c(i)
# [maximise], where dist_E(i) = BFS distance from the end node E and
# c(i) = count of PREACTIVE+ACTIVE neighbours.
#
# Weights are auto-scaled per component:
#   eff_w1 = w1 * max_degree,  eff_w2 = w2 * max_dist_E
# so both terms span comparable ranges regardless of mesh order.
#
# Status: INACTIVE=0, PREACTIVE=1, ACTIVE=2, NUMBERED=3.
# ===================================================================
function _sloan_perm(adj::Table; w1::Int=1, w2::Int=2)
  nfree   = length(adj)
  degrees = diff(adj.ptrs)

  queue = Vector{Int}(undef, nfree)
  level = Vector{Int}(undef, nfree)
  state = zeros(Int, nfree)
  gen   = Ref(0)

  function bfs!(start)
    cg = (gen[] += 1)
    state[start] = cg; level[start] = 0
    queue[1] = start; qhead = qtail = 1
    max_level = 0
    while qhead <= qtail
      v = queue[qhead]; qhead += 1; lv = level[v]
      for k in datarange(adj, v)
        u = adj.data[k]
        if state[u] % UInt < cg % UInt
          state[u] = cg; level[u] = lv + 1
          lv + 1 > max_level && (max_level = lv + 1)
          qtail += 1; queue[qtail] = u
        end
      end
    end
    last = qtail
    while last > 1 && level[queue[last-1]] == max_level; last -= 1; end
    return qtail, last
  end

  dist_E  = Vector{Int}(undef, nfree)
  f       = fill(typemin(Int), nfree)
  c       = zeros(Int, nfree)
  status  = zeros(Int8, nfree)

  heap = Tuple{Int,Int}[]

  function pq_push!(node, pr)
    f[node] = pr
    push!(heap, (pr, node))
    i = length(heap)
    @inbounds while i > 1
      p = i >> 1; heap[p][1] >= heap[i][1] && break
      heap[p], heap[i] = heap[i], heap[p]; i = p
    end
  end

  function pq_siftdown!()
    i = 1; n = length(heap)
    @inbounds while true
      c_idx = 2i
      c_idx > n && break
      c_idx + 1 <= n && heap[c_idx+1][1] > heap[c_idx][1] && (c_idx += 1)
      heap[c_idx][1] <= heap[i][1] && break
      heap[c_idx], heap[i] = heap[i], heap[c_idx]; i = c_idx
    end
  end

  function pq_pop!()
    @inbounds while !isempty(heap)
      pr, node = heap[1]
      if length(heap) == 1; pop!(heap)
      else heap[1] = pop!(heap); pq_siftdown!()
      end
      f[node] == pr && return node
    end
    return 0
  end

  perm = Vector{Int}(undef, nfree)
  idx  = 1

  while idx <= nfree
    seed, seed_deg = 0, typemax(Int)
    for i in 1:nfree
      status[i] == 0 && degrees[i] < seed_deg && (seed_deg = degrees[i]; seed = i)
    end

    E = _george_liu_peripheral!(bfs!, queue, level, degrees, seed)
    n_comp, last_e = bfs!(E)

    S, S_deg = 0, typemax(Int)
    max_dist = 0; max_deg = 0
    for qi in 1:n_comp
      v = queue[qi]
      dist_E[v] = level[v]
      level[v] > max_dist && (max_dist = level[v])
      degrees[v] > max_deg && (max_deg = degrees[v])
      qi >= last_e && degrees[v] < S_deg && (S_deg = degrees[v]; S = v)
    end
    eff_w1 = w1 * max(max_deg, 1)
    eff_w2 = w2 * max(max_dist, 1)

    status[S] = 1; c[S] = 0; pq_push!(S, eff_w1 * dist_E[S])

    comp_done = 0
    while comp_done < n_comp
      v = pq_pop!(); v == 0 && break

      if status[v] == 1
        status[v] = 2
        for k in datarange(adj, v)
          u = adj.data[k]
          if status[u] == 0
            status[u] = 1
            c_u = 0
            for k2 in datarange(adj, u)
              w = adj.data[k2]
              if status[w] == 1 || status[w] == 2
                c_u += 1
                c[w] += 1
                pq_push!(w, eff_w1 * dist_E[w] + eff_w2 * c[w])
              end
            end
            c[u] = c_u
            pq_push!(u, eff_w1 * dist_E[u] + eff_w2 * c_u)
          end
        end
        pq_push!(v, eff_w1 * dist_E[v] + eff_w2 * c[v])

      elseif status[v] == 2
        status[v] = 3; f[v] = typemin(Int)
        perm[idx] = v; idx += 1; comp_done += 1
      end
    end

    for qi in 1:n_comp; state[queue[qi]] = -1; end
  end

  iperm = Vector{Int}(undef, nfree)
  for i in 1:nfree; iperm[perm[i]] = i; end
  return iperm
end

# ===================================================================
# Coordinate-based permutation (Lagrangian spaces).
# Sorts free DOFs by their spatial coordinates.
# by: key function applied to each VectorValue before comparison.
#     Defaults to Tuple, giving lexicographic (x,y,z) ordering.
# lt: less-than predicate used by sortperm (default isless).
# ===================================================================
function _coordinates_perm(space::FESpace; by=Tuple, lt=isless)
  coords = get_free_dof_coordinates(space)
  invperm(sortperm(coords; by=by, lt=lt))
end

# ===================================================================
# Main public API
# ===================================================================

"""
    compute_dof_permutation(space::FESpace, algorithm::Symbol=:rcm; adj, by, lt, kwargs...)

Compute a free-DOF permutation `iperm` where `iperm[old_dof] = new_dof`.

`algorithm` selects the reordering strategy:
- `:rcm` — Reverse Cuthill-McKee with George-Liu pseudo-peripheral node
  selection. Best for bandwidth and profile reduction on structured meshes.
- `:sloan` — Sloan's algorithm. Minimises the active wavefront; may
  outperform `:rcm` on profile for irregular or 3D meshes.
  Accepts keyword arguments `w1::Int=1` and `w2::Int=2` to tune the
  distance/connectivity weight balance.
- `:coordinates` — Sort free DOFs by their spatial coordinates (Lagrangian
  spaces only). Keyword `by` maps each `VectorValue` to a sortable key
  (default `Tuple`, giving lexicographic x-y-z ordering); `lt` is the
  less-than predicate (default `isless`). Example: `by=x->x[1]` sorts by
  x-coordinate only.

For `:rcm` and `:sloan`, the optional keyword argument `adj` is the DOF
adjacency `Table` returned by `compute_adjacency`. It defaults to the graph
built from the cell-DOF connectivity of `space`. Pass a pre-built table to
avoid redundant work when calling multiple algorithms on the same space.
Passing `adj` when `algorithm == :coordinates` has no effect.
"""
function compute_dof_permutation(space::FESpace, algorithm::Symbol=:rcm; adj=nothing, by=Tuple, lt=isless, kwargs...)
  if algorithm == :coordinates
    return _coordinates_perm(space; by=by, lt=lt)
  end
  if isnothing(adj)
    adj = compute_adjacency(get_cell_dof_ids(space), num_free_dofs(space))
  end
  @check length(adj) == num_free_dofs(space)
  algorithm == :rcm   && return _rcm_perm(adj)
  algorithm == :sloan && return _sloan_perm(adj; kwargs...)
  error("Unknown algorithm :$algorithm. Supported: :rcm, :sloan, :coordinates")
end

"""
    reindex_free_dof_ids(space::FESpace, algorithm::Symbol=:rcm; kwargs...)

Return a new `FESpace` with free DOF IDs reordered using `algorithm`.
See [`compute_dof_permutation`](@ref) for supported algorithms and keyword arguments.
Dirichlet DOF IDs are left unchanged.

# Examples
```julia
model = CartesianDiscreteModel((0,1,0,1), (8,8))
V = FESpace(model, ReferenceFE(lagrangian, Float64, 2); dirichlet_tags="boundary")

V_rcm    = reindex_free_dof_ids(V, :rcm)
V_coords = reindex_free_dof_ids(V, :coordinates)
V_xonly  = reindex_free_dof_ids(V, :coordinates; by=x->x[1])
```
"""
function reindex_free_dof_ids(space::FESpace, algorithm::Symbol=:rcm; kwargs...)
  reindex_free_dof_ids(space, compute_dof_permutation(space, algorithm; kwargs...))
end
