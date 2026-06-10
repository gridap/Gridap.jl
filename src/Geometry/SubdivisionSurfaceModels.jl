# Loop surface model

function loop_surface_model(chart_model::DiscreteModel{2}, geo_map)
  @check all(==(TRI), get_polytopes(chart_model)) "Loop surface require surface mesh of triangles."
  @check all(get_faces(get_grid_topology(chart_model),0,1) .|> length .|> ==(6) ) """
    Loop surface is only implemented on boundary-less meshes with regular vertices (valence 6).
  """
  @check isa(get_grid_topology(chart_model), UnstructuredGridTopology)

  # Compute the correct physical node

  geo_map_field = GenericField(geo_map)
  chart_node_coords = get_node_coordinates(chart_model)
  node_coords = evaluate(geo_map_field, chart_node_coords)

  loop_surface_model(chart_model::DiscreteModel{2}, node_coords)
end

function loop_surface_model(chart_model::DiscreteModel{2}, node_coordinates::AbstractVector{<:Point})
  @check all(==(TRI), get_polytopes(chart_model)) "Loop surface require surface mesh of triangles."
  @check all(get_faces(get_grid_topology(chart_model),0,1) .|> length .|> ==(6) ) """
    Loop surface is only implemented on boundary-less meshes with regular vertices (valence 6).
  """
  @check isa(get_grid_topology(chart_model), UnstructuredGridTopology)

  # Compute physical vertices
  vertex_to_node = get_vertex_node(chart_model)
  vertex_coords = evaluate(Reindex(node_coordinates), vertex_to_node)
  Dp = length(first(vertex_coords))
  T = eltype(first(vertex_coords))

  # Compute the cell_map

  ptopo_cell = PatchTopology(chart_model)
  loop_ptopo = extend_patches_by_single_layer(ptopo_cell)
  n_patchs = num_patches(loop_ptopo)
  @check n_patchs == num_cells(chart_model)

  TRI_vertex_coords = get_vertex_coordinates(TRI)
  generating_splines = ReferenceFEs._box_splines_222(T, TRI_vertex_coords)
  p_gen_splines = Fill(generating_splines , n_patchs)
  p_vertices = get_patch_faces(loop_ptopo, 0)
  p_coords = lazy_map(Broadcasting(Reindex(vertex_coords)), p_vertices)
  cell_map = lazy_map(linear_combination, p_coords, p_gen_splines)

  # same grid, topology and labeling with new cell_map and vertex/nodes coordinates

  cgrid = get_grid(chart_model)
  grid = UnstructuredGrid(
    node_coordinates,
    get_cell_node_ids(cgrid),
    get_reffes(cgrid), # TODO they make no sense
    get_cell_type(cgrid),
    OrientationStyle(cgrid),
    nothing, # get_facet_normal(cgrid) # TODO compute them if they were in chart_model?
    cell_map
  )

  ctopo = get_grid_topology(chart_model)
  orien = OrientationStyle(ctopo)
  topo = UnstructuredGridTopology{2,Dp,T,typeof(orien)}(
    vertex_coords,
    ctopo.n_m_to_nface_to_mfaces,
    ctopo.cell_type,
    ctopo.polytopes,
    orien
  )

  UnstructuredDiscreteModel(grid, topo, chart_model.face_labeling)
end

