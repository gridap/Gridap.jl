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

`topo` is the topology of the control grid. It must only have isolated irregular
vertices, that is, if a vertex has valence ≠6, all it's neighboors at topological
distance up to 2 are of valence 6.

The map is evaluated at the vertex vector of a triangle, and return the vertices
of it's Loop patch in the following expected order:

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
  loop_vertices, v_to_av_caches = cache

  setsize!(loop_vertices, (12,))
  compute_loop_patch_vertices!(
    loop_vertices,  element_vertices, k.vertex_to_adj_vertices, v_to_av_caches
  )

  loop_vertices
end

function compute_loop_patch_vertices!(
  patch_vertices, q_vertices, vertex_to_vertices, v_to_av_caches
)

  v_to_av_cache1, v_to_av_cache2 = v_to_av_caches
  Base.@propagate_inbounds adjascent_vertices(vi) = getindex!(v_to_av_cache2, vertex_to_vertices, vi)
  Base.@propagate_inbounds function find_adjto_adjto_butnot(v1, v2, v3)
    v1_adj = getindex!(v_to_av_cache1, vertex_to_vertices, v1)
    findfirstvalue(v1_adj) do vi
      vi ≠ v3  &&
      vi ∈ adjascent_vertices(v2)
    end
  end

  # element vertices
  #v4, v4iq = findmin(q_vertices)
  #v8, v8iq = findmax(q_vertices)
  #v7iq = findfirst(!in((v4iq,v8iq)), 1:3)
  #v7 = q_vertices[v7iq]
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

