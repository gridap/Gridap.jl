# Loop surface model

"""
    loop_surface_model(chart_model::DiscreteModel{2}, geo_map)

Loop surface model constructor. Uses a base `chart_model`, and a given
geometrical map

    `geo_map` : `chart_model` -> S,

and return a model approximating S using a Loop surface such that:
- the control mesh topology and labels matches that of `chart_model`,
- the control vertices coordinates are found using a L2 projection of `geo_map` onto the Loop spline space.
"""
function loop_surface_model(chart_model::DiscreteModel{2}, geo_map::Union{Function, Field})
  @check all(==(TRI), get_polytopes(chart_model)) "Loop surface require surface mesh of triangles."
  @check all(get_faces(get_grid_topology(chart_model),0,1) .|> length .|> ==(6) ) """
    Loop surface is only implemented on boundary-less meshes with regular vertices (valence 6).
  """
  @check isa(get_grid_topology(chart_model), UnstructuredGridTopology)

  # get (loosely) approximate control vertex coordinates by interpolating the geo_map
  geo_map_field = GenericField(geo_map)
  chart_vertex_coordinates = get_vertex_coordinates(get_grid_topology(chart_model))
  geo_map_cache = return_cache(geo_map_field, chart_vertex_coordinates)
  vertex_coordinates = evaluate!(geo_map_cache, geo_map_field, chart_vertex_coordinates)

  # Definition of a Loop FE space
  P = eltype(vertex_coordinates)
  loop_reffe = LoopRefFE(P,TRI)
  Ω_chart = Triangulation(chart_model)
  dΩ = Measure(Ω_chart, 8)
  # this dispatches to hacky subdivision_surfaces/src/FESpaces/LoopFESpaces.jl
  U = FESpace(Ω_chart, loop_reffe)

  # L2 projection
  a(Φ, Ψ) = ∫(Ψ ⋅ Φ)dΩ
  l(   Ψ) = ∫(Ψ ⋅ geo_map)dΩ
  op = AffineFEOperator(a, l, U, U)
  Φh = solve(op)
  e = Φh-geo_map
  errl2 = sqrt(sum(a(e,e)))

  @info "Loop surface L²-projection L²-error: $errl2"

  # Retrieving the control vertices coordinates from the DoF vector
  D = length(P)
  n_vertices = length(vertex_coordinates)
  M = reshape(Φh.free_values, (D, n_vertices));
  fitted_control_vertices = collect( P(r...) for r in eachcol(M));

  loop_surface_model(get_grid_topology(chart_model), fitted_control_vertices)
end

"""
    loop_surface_model(
      control_topo::GridTopology{2},
      control_vertex_coordinates::AbstractVector{<:Point}),
      face_labeling = FaceLabeling(topo)
    )

Lower level Loop surface model constructor, taking a control mesh
topology, the associated control vertices coordinates, and optional face labeling.
"""
function loop_surface_model(
  topo::GridTopology{2},
  control_vertex_coordinates::AbstractVector{<:Point},
  face_labeling = FaceLabeling(topo)
)
  @check all(==(TRI), get_polytopes(topo)) "Loop surface require surface mesh of triangles."
  @check all(get_faces(topo,0,1) .|> length .|> ==(6) ) """
    Loop surface is only implemented on boundary-less meshes with regular vertices (valence 6).
  """
  @check isa(topo, UnstructuredGridTopology)
  Dp = length(first(control_vertex_coordinates))
  T = eltype(first(control_vertex_coordinates))

  # cell_map computation

  # Define the global control vertices coordinates of each Loop patch (for each element/triangle)
  n_patchs = num_cells(topo)
  cell_to_vertices = get_faces(topo, 2, 0)
  patch_loop_vertices = Table(lazy_map(LoopPatchVerticesMap(topo), cell_to_vertices))
  #println.(patch_loop_vertices)
  p_coords = lazy_map(Broadcasting(Reindex(control_vertex_coordinates)), patch_loop_vertices)

  # Create the Loop generating splines basis, and the cell maps by compination with the vertex coordinates
  TRI_vertex_coords = get_vertex_coordinates(TRI)
  generating_splines = _box_splines_222(T, TRI_vertex_coords)
  p_gen_splines = Fill(generating_splines, n_patchs)
  cell_map = lazy_map(linear_combination, p_coords, p_gen_splines)

  # grid with trivial reffes topology and labeling with new cell_map and vertex/nodes coordinates
  reffes = [ LagrangianRefFE(Point{Dp,T}, TRI, 1) ]
  orien = OrientationStyle(topo)
  cell_types = collect(Fill(Int8(1), num_cells(topo)))
  loop_grid = UnstructuredGrid(
    control_vertex_coordinates,
    cell_to_vertices, reffes, cell_types,
    orien, nothing, cell_map # TODO compute normals ?
  )

  loop_topo = UnstructuredGridTopology{2,Dp,T,typeof(orien)}(
    control_vertex_coordinates,
    topo.n_m_to_nface_to_mfaces,
    topo.cell_type,
    topo.polytopes,
    orien
  )

  UnstructuredDiscreteModel(loop_grid, loop_topo, face_labeling)
end

