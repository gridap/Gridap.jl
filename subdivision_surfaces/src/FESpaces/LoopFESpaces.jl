# FE Space whose elements's polytope faces do not own all the DOFs of it's ref-FE

function FESpace(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{<:Loop}};
  conformity=nothing,
  trian = Triangulation(model),
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type=nothing,
  scale_dof=false,
  global_meshsize=nothing
)

  conf = Conformity(testitem(cell_reffe),conformity)
  cell_fe = FESpaces.CellFE(model, cell_reffe, conf; scale_dof, global_meshsize)
  fespace = FESpace(
    model, cell_fe; trian, labels, dirichlet_tags, dirichlet_masks, constraint, vector_type
  )

  #@assert fespace.

  reffe = testitem(cell_reffe)
  shapefuns = get_shapefuns(reffe)
  L = num_indep_components(value_type(shapefuns))

  ctopo = get_grid_topology(model)
  patch_vertex_map = LoopPatchVerticesMap(ctopo)
  cell_to_vertices = get_cell_vertices(ctopo)
  cell_patch_vertices = Table(lazy_map(patch_vertex_map, cell_to_vertices))

  for k in eachindex(cell_patch_vertices)
    k_dofs_ids = @view fespace.cell_dofs_ids[k]
    k_loop_ids = @view cell_patch_vertices[k]
    k_dof = 1
    for v in k_loop_ids
      for ll in 1:L
        k_dofs_ids[k_dof] = L*(v-1) + ll
        k_dof += 1
      end
    end
  end

  fespace
end


########################
# LoopPatchVerticesMap #
########################

"""
    LoopPatchVerticesMap(topo::GridTopology)

`topo` is the topology of the control grid. It must only have isolated extraordinary
vertices (EVs), that is, if a vertex has valence ≠6, all it's neighboors at
topological distance up to 2 are of valence 6.

The map is evaluated at the vertex vector of a triangle, and returns the
vertices of it's Loop patch. See [`compute_loop_patch_vertices!`](@ref) for the orderings.
"""
struct LoopPatchVerticesMap <: Map
  vertex_to_adj_vertices::Table{Int32,Vector{Int32},Vector{Int32}}
end

function LoopPatchVerticesMap(topo::GridTopology)
  node_adj_graph = compute_graph(topo, 0, 1; self_loops=false)
  vertex_to_adj_vertices = Table(node_adj_graph.rowval, node_adj_graph.colptr)
  LoopPatchVerticesMap(vertex_to_adj_vertices)
end

function return_cache(k::LoopPatchVerticesMap, element_vertices)
  # I need two caches because I lookup for another one while iterating over the
  # first one in find_adjto_adjto_butnot below, so the outer one shouln't be mutated
  v_to_av_cache1 = array_cache(k.vertex_to_adj_vertices)
  v_to_av_cache2 = array_cache(k.vertex_to_adj_vertices)
  loop_vertices_cache = CachedArray(Int32[])
  loop_vertices_cache, (v_to_av_cache1, v_to_av_cache2)
end

function evaluate!(cache, k::LoopPatchVerticesMap, element_vertices)
  loop_vertices, (v_to_av_cache1, v_to_av_cache2) = cache

  Base.@propagate_inbounds adjascent_vertices(vi) = getindex!(v_to_av_cache2, k.vertex_to_adj_vertices, vi)
  # Finds a vertex adjacent to v1 and to v2, other than v3, and, if given,
  # not adjacent to v4 either.
  Base.@propagate_inbounds function find_adjto_adjto_butnot(v1, v2, v3, v4=nothing)
    v1_adj = getindex!(v_to_av_cache1, k.vertex_to_adj_vertices, v1)
    findfirstvalue(v1_adj) do vi
      vi ≠ v3  &&
      vi ∈ adjascent_vertices(v2) &&
      (isnothing(v4) || vi ∉ adjascent_vertices(v4))
    end
  end

  # local id (1, 2 or 3) of the extraordinary vertex (EV) within
  # element_vertices, if any; defaulted to 1 (unused) for a regular triangle
  l_EV_id = findfirst(v -> length(adjascent_vertices(v)) ≠ 6, element_vertices)
  N = isnothing(l_EV_id) ? 6 : length(adjascent_vertices(element_vertices[l_EV_id]))

  setsize!(loop_vertices, (N+6,))
  compute_loop_patch_vertices!(
    loop_vertices, element_vertices, something(l_EV_id, 1), Val(N),
    adjascent_vertices, find_adjto_adjto_butnot
  )

  loop_vertices
end

"""
    compute_loop_patch_vertices!(patch_vertices, q_vertices, l_EV_id, ::Val{6}, ...)

Fills `patch_vertices` (length 12) with the global vertices corresponding to the
following regular Loop patch local vertex ids. 4-7-8 correspond to the vertices
`q_vertices` (`l_EV_id` unused)

```plain
       (1)-------(2)
      /   \\     /  \\
     /     \\   /    \\
    (3)-----(4)------(5)
   /   \\   /   \\   /   \\
  /     \\ /     \\ /     \\
(6)-----(7)-----(8)-----(9)
  \\    /  \\    /  \\    /
   \\  /    \\  /    \\  /
   (10)----(11)----(12)
```
"""
function compute_loop_patch_vertices!(
  patch_vertices, q_vertices, l_EV_id, ::Val{6},
  adjascent_vertices, find_adjto_adjto_butnot
)

  v4, v7, v8 = q_vertices
  msg =  "only regular patch are implemented (vertices have valence 6), no boundary"
  @check 6 == length(adjascent_vertices(v4)) msg
  @check 6 == length(adjascent_vertices(v7)) msg
  @check 6 == length(adjascent_vertices(v8)) msg

  @inbounds begin
    # vertices adjascent to two element vertices
    v3  = find_adjto_adjto_butnot(v4, v7, v8)
    v11 = find_adjto_adjto_butnot(v7, v8, v4)
    v5  = find_adjto_adjto_butnot(v8, v4, v7)

    # remaining vertices
    patch_vertices[1]  = find_adjto_adjto_butnot(v3, v4, v7)
    patch_vertices[2]  = find_adjto_adjto_butnot(v4, v5, v8)
    patch_vertices[3]  = v3
    patch_vertices[4]  = v4
    patch_vertices[5]  = v5
    patch_vertices[6]  = find_adjto_adjto_butnot(v3, v7, v4)
    patch_vertices[7]  = v7
    patch_vertices[8]  = v8
    patch_vertices[9]  = find_adjto_adjto_butnot(v5, v8, v4)
    patch_vertices[10] = find_adjto_adjto_butnot(v7, v11,v8)
    patch_vertices[11] = v11
    patch_vertices[12] = find_adjto_adjto_butnot(v11,v8, v7)

    @check all(patch_vertices[ [1 2 3 5 6 9 10 11 12] ] .|> adjascent_vertices .|> length .|> ==(6)) msg
  end

  patch_vertices
end

"""
    compute_loop_patch_vertices!(patch_vertices, q_vertices, l_EV_id, ::Val{N}, ...) where N

Fills `patch_vertices` (length `N+6`) with the vertices of the irregular Loop
patch of triangle `q_vertices`, whose `l_EV_id`-th vertex is the extraordinary
vertex (EV) of valence `N`, in the numbering of Stam (1998), Fig. 2:
- `1` is the EV, `2` and `N+1` are the two other vertices of `q`
- `2,…,N+1` are its `N` neighbors (`7,…,N-1` are elided since they depend on `N`)
- `N+2,…,N+6` complete the 1-neighborhood patch.

Example of patches:
```plain
N > 6:                           N = 3:

5------6 ···· N------N+6                          N
|\\     |     /|     /|                         / /|\\ \\
| \\    |   /  |   /  |                      /    /|\\    \\
|   \\  |  /   |  /   |                   /      / | \\      \\
|     \\|/     |/     |                /        /  |  \\        \\
4------1------N+1----N+5           N+4        /   1   \\          N+6
|     /|     /|     /             /  \\        /  / \\  \\        /  \\
|   /  |   /  |   /               /   \\      / /     \\ \\      /   \\
|  /   |  /   |  /                /     \\   / /       \\ \\   /     \\
|/     |/     |/                 /       \\  /           \\  /       \\
3------2------N+2               /          2-------------N+1        \\
|     /|     /                  /       /   \\           /   \\       \\
|   /  |   /                    /     /       \\       /       \\     \\
|  /   |  /                    /   /           \\     /           \\   \\
|/     |/                     /  /               \\ /               \\  \\
N+4----N+3                    N+3-----------------N+2-----------------N+5
```
"""
function compute_loop_patch_vertices!(
  patch_vertices, q_vertices, l_EV_id, ::Val{N},
  adjascent_vertices, find_adjto_adjto_butnot
) where N

  @inbounds begin
    ev = q_vertices[l_EV_id]
    v2 = q_vertices[mod1(l_EV_id+1, 3)]
    vN1 = q_vertices[mod1(l_EV_id+2, 3)]

    patch_vertices[1]   = ev
    patch_vertices[2]   = v2
    patch_vertices[N+1] = vN1

    msg = """
      irregular Loop patch requires an isolated extraordinary vertex: the two
      other triangle vertices, and the whole ring around the extraordinary one,
      must be regular (v2lence 6), no boundary.
    """
    @check N == length(adjascent_vertices(ev)) msg
    @check N ≠ 6 msg
    @check 6 == length(adjascent_vertices(v2)) msg
    @check 6 == length(adjascent_vertices(vN1)) msg

    # walk the ring from 2 to N+1 to fill them
    prev2, prev1 = vN1, v2
    for i in 3:N
      vi = find_adjto_adjto_butnot(ev, prev1, prev2)
      patch_vertices[i] = vi
      prev2, prev1 = prev1, vi
    end

    v3 = patch_vertices[3]
    vN = patch_vertices[N]
    vN2 = find_adjto_adjto_butnot(v2, vN1, ev, ev)

    patch_vertices[N+2] = vN2
    patch_vertices[N+3] = find_adjto_adjto_butnot(v2,  vN2, vN1, ev)
    patch_vertices[N+4] = find_adjto_adjto_butnot(v3,  v2,  ev,  ev)
    patch_vertices[N+5] = find_adjto_adjto_butnot(vN2, vN1, v2,  ev)
    patch_vertices[N+6] = find_adjto_adjto_butnot(vN,  vN1, ev,  ev)

    for i in 3:N+6
      @check 6 == length(adjascent_vertices(patch_vertices[i])) msg
    end
  end

  patch_vertices
end
