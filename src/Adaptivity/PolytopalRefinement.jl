
# Notes: 
# We need to restrict the input to PolytopalDiscreteModel. This would not be strictly 
# necessary (we can still aggregate other models), but we are currently limited 
# because we need that all polytopes are positively-oriented rotation systems (like GeneralPolytopes). 
# Unfortunately our ExtrusionPolytopes are not.
# To make it work for any model, we would have to reindex the face nodes to ensure the faces
# are correctly oriented.

function coarsen(model::Geometry.PolytopalDiscreteModel,ptopo::Geometry.PatchTopology)
  new_polys, new_connectivity = generate_patch_polytopes(model,ptopo)

  vertex_coords = get_vertex_coordinates(get_grid_topology(model))
  new_to_old = unique(new_connectivity.data)
  old_to_new = find_inverse_index_map(new_to_old)
  new_vertex_coords = vertex_coords[new_to_old]
  map!(old -> old_to_new[old], new_connectivity.data, new_connectivity.data)

  new_topo = Geometry.PolytopalGridTopology(new_vertex_coords,new_connectivity,new_polys)
  new_grid = Geometry.PolytopalGrid(new_topo)
  new_labels = FaceLabeling(new_topo)
  new_model = Geometry.PolytopalDiscreteModel(new_grid,new_topo,new_labels)

  return new_model
end

function generate_patch_polytopes(model::DiscreteModel, ptopo::Geometry.PatchTopology)
  topo = get_grid_topology(model)
  @assert topo === ptopo.topo

  D = num_cell_dims(topo)
  polys = get_polytopes(topo)
  vertex_coordinates = get_vertex_coordinates(topo)
  cell_to_vertices = get_faces(topo,D,0)
  face_to_vertices = get_faces(topo,D-1,0)

  npatches = Geometry.num_patches(ptopo)
  patch_cells = Geometry.get_patch_cells(ptopo)
  patch_boundary = Geometry.PatchBoundaryTriangulation(model,ptopo)
  tface_to_face = patch_boundary.trian.trian.tface_to_mface
  patch_to_tfaces = patch_boundary.glue.patch_to_tfaces
  face_glue = patch_boundary.trian.glue

  new_connectivity = Vector{Vector{Int32}}(undef, npatches)
  new_polys = Vector{GeneralPolytope{D}}(undef, npatches)
  Threads.@threads for patch in 1:npatches
    cells = view(patch_cells,patch)
    if isone(length(cells))
      new_polys[patch] = polys[first(cells)]
      new_connectivity[patch] = cell_to_vertices[first(cells)]
      continue
    end
    
    _vertices = Set{Int32}()
    tfaces = view(patch_to_tfaces,patch)
    connectivity = Vector{Vector{Int32}}(undef, length(tfaces))
    for (k,tface) in enumerate(tfaces)
      face = tface_to_face[tface]
      cell = face_glue.face_to_cell[tface]
      lface = face_glue.face_to_lface[tface]
      pindex = face_glue.cell_to_lface_to_pindex[cell][lface]
      fperm = get_vertex_permutations(Polytope{D-1}(polys[cell],lface))[pindex]

      connectivity[k] = face_to_vertices[face][fperm]
      push!(_vertices, connectivity[k]...)
    end

    vertices = collect(_vertices)
    coords = vertex_coordinates[vertices]
    glob_to_loc = Dict{Int32,Int32}( v => k for (k,v) in enumerate(vertices))
    foreach(y -> map!(x -> glob_to_loc[x], y, y), connectivity)

    new_polys[patch], vperm = ReferenceFEs.polytope_from_faces(D,coords,connectivity)
    new_connectivity[patch] = vertices[vperm]
  end

  return new_polys, Table(new_connectivity)
end
