# Loop surface model

function loop_surface_model(chart_model::DiscreteModel{2}, geo_map)
  @check all(==(TRI), get_polytopes(chart_model)) "Loop surface require surface mesh of triangles."
  @check all(get_faces(get_grid_topology(chart_model),0,1) .|> length .|> ==(6) ) """
    Loop surface is only implemented on boundary-less meshes with regular vertices (valence 6).
  """
  @check isa(get_grid_topology(chart_model), UnstructuredGridTopology)

  println("Computing approximate control node")

  geo_map_field = GenericField(geo_map)
  chart_vertex_coordinates = get_vertex_coordinates(get_grid_topology(chart_model))
  vertex_coordinates = evaluate(geo_map_field, chart_vertex_coordinates )

  loop_surface_model(chart_model::DiscreteModel{2}, vertex_coordinates)
end

"""
    loop_surface_model(chart_model::DiscreteModel{2}, geo_map)
    loop_surface_model(chart_model::DiscreteModel{2},
      control_node_coordinates::AbstractVector{<:Point})


Builds a model for a Loop subdivision surface
"""
function loop_surface_model(
  chart_model::DiscreteModel{2},
  control_vertex_coordinates::AbstractVector{<:Point}
  #control_node_coordinates::AbstractVector{<:Point}
)
  @check all(==(TRI), get_polytopes(chart_model)) "Loop surface require surface mesh of triangles."
  @check all(get_faces(get_grid_topology(chart_model),0,1) .|> length .|> ==(6) ) """
    Loop surface is only implemented on boundary-less meshes with regular vertices (valence 6).
  """
  @check isa(get_grid_topology(chart_model), UnstructuredGridTopology)
  Dp = length(first(control_vertex_coordinates))
  T = eltype(first(control_vertex_coordinates))

  # physical/control vertices computation
  #vertex_to_node = get_vertex_node(chart_model)
  #control_vertex_coords = evaluate(Reindex(control_node_coordinates), vertex_to_node)

  # cell_map computation

  # Define the control vertices coordinates of each Loop patch (for each element/triangle)
  n_patchs = num_cells(chart_model)
  ctopo = get_grid_topology(chart_model)
  cell_to_vertices = get_faces(ctopo, 2, 0)
  patch_loop_vertices = Table(lazy_map(LoopPatchVerticesMap(ctopo), cell_to_vertices))
  #println.(patch_loop_vertices)
  p_coords = lazy_map(Broadcasting(Reindex(control_vertex_coordinates)), patch_loop_vertices)

  # Create the Loop generating splines basis, and the cell maps by compination with the vertex coordinates
  TRI_vertex_coords = get_vertex_coordinates(TRI)
  generating_splines = _box_splines_222(T, TRI_vertex_coords)
  p_gen_splines = Fill(generating_splines, n_patchs)
  cell_map = lazy_map(linear_combination, p_coords, p_gen_splines)

  # same grid, topology and labeling with new cell_map and vertex/nodes coordinates
  cgrid = get_grid(chart_model)
  #reffes = [LoopRefFE(Point{Dp,T}, TRI)]
  grid = UnstructuredGrid(
    #control_node_coordinates,
    control_vertex_coordinates,
    #get_cell_node_ids(cgrid),
    cell_to_vertices,
    get_reffes(cgrid),
    #reffes,
    get_cell_type(cgrid),
    OrientationStyle(cgrid),
    nothing, # get_facet_normal(cgrid) # TODO compute them if they were in chart_model?
    cell_map
  )

  ctopo = get_grid_topology(chart_model)
  orien = OrientationStyle(ctopo)
  topo = UnstructuredGridTopology{2,Dp,T,typeof(orien)}(
    control_vertex_coordinates,
    ctopo.n_m_to_nface_to_mfaces,
    ctopo.cell_type,
    ctopo.polytopes,
    orien
  )

  #GenericDiscreteModel(grid, topo, chart_model.face_labeling)
  UnstructuredDiscreteModel(grid, topo, chart_model.face_labeling)
end


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

  v_to_av_cache1,  v_to_av_cache2 = v_to_av_caches
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

