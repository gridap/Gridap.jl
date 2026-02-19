
# Notes: 
# We need to restrict the input to PolytopalDiscreteModel. This would not be strictly 
# necessary (we can still aggregate other models), but we are currently limited 
# because we need that all polytopes are positively-oriented rotation systems (like GeneralPolytopes). 
# Unfortunately our ExtrusionPolytopes are not.
# To make it work for any model, we would have to reindex the face nodes to ensure the faces
# are correctly oriented.

function coarsen(model::Geometry.PolytopalDiscreteModel,ptopo::Geometry.PatchTopology; return_glue=false)
  @check Geometry.is_partition(ptopo) "The patch topology is not a valid partition of the model"
  new_polys, new_connectivity = generate_patch_polytopes(model,ptopo)

  vertex_coords = get_vertex_coordinates(get_grid_topology(model))
  new_to_old_nodes = unique(new_connectivity.data)
  old_to_new_nodes = find_inverse_index_map(new_to_old_nodes)
  new_vertex_coords = vertex_coords[new_to_old_nodes]
  map!(old -> old_to_new_nodes[old], new_connectivity.data, new_connectivity.data)

  new_topo = Geometry.PolytopalGridTopology(new_vertex_coords,new_connectivity,new_polys)
  new_grid = Geometry.PolytopalGrid(new_topo)
  new_labels = FaceLabeling(new_topo)
  new_model = Geometry.PolytopalDiscreteModel(new_grid,new_topo,new_labels)

  if !(return_glue)
    return new_model
  else
    glue = generate_patch_adaptivity_glue(
      ptopo, get_grid_topology(model), new_topo, old_to_new_nodes, new_to_old_nodes,
    )
    return new_model, glue
  end
end

function generate_patch_adaptivity_glue(
  ptopo, ftopo, ctopo, fine_to_coarse_nodes, coarse_to_fine_nodes
)
  Dc = num_cell_dims(ftopo)
  coarse_to_fine_cells = Geometry.get_patch_cells(ptopo)
  fine_to_coarse_cells = Arrays.flatten_partition(coarse_to_fine_cells, num_faces(ftopo,Dc))
  generate_patch_adaptivity_glue(
    ftopo, ctopo, fine_to_coarse_cells, coarse_to_fine_cells, fine_to_coarse_nodes, coarse_to_fine_nodes
  )
end

function generate_patch_adaptivity_glue(
  ftopo, ctopo, fine_to_coarse_cells, coarse_to_fine_cells, fine_to_coarse_nodes, coarse_to_fine_nodes
)
  Dc = num_cell_dims(ftopo)
  fine_to_coarse_faces = Vector{Vector{Int32}}(undef, Dc+1)
  fine_to_coarse_faces[1] = fine_to_coarse_nodes
  fine_to_coarse_faces[Dc+1] = fine_to_coarse_cells
  for d in 1:Dc-1
    cface_to_cnodes = Geometry.get_faces(ctopo,d,0)
    fnode_to_ffaces = Geometry.get_faces(ftopo,0,d)
    coarse_to_fine_faces = zeros(Int32, num_faces(ctopo,d))
    for cface in eachindex(cface_to_cnodes)
      cnodes = view(cface_to_cnodes, cface)
      fnodes = view(coarse_to_fine_nodes, cnodes)
      if all(!iszero, fnodes)
        fface = only(intersect((view(fnode_to_ffaces,fnode) for fnode in fnodes)...))
        coarse_to_fine_faces[cface] = fface
      end
    end
    fine_to_coarse_faces[d+1] = Arrays.find_inverse_index_map(coarse_to_fine_faces, num_faces(ftopo,d))
  end

  fine_child_ids = Arrays.find_local_index(fine_to_coarse_cells, coarse_to_fine_cells)
  refinement_rules = Fill(WhiteRefinementRule(TRI), length(coarse_to_fine_cells))
  is_refined = select_refined_cells(fine_to_coarse_cells)
  return AdaptivityGlue(
    RefinementGlue(), fine_to_coarse_faces, fine_child_ids, refinement_rules, is_refined, coarse_to_fine_cells
  )
end

function generate_patch_polytopes(
  model::DiscreteModel{D}, ptopo::Geometry.PatchTopology{D}
) where D
  topo = get_grid_topology(model)
  @assert topo === ptopo.topo

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
    
    nv = 0
    glob_to_loc = OrderedDict{Int32,Int32}()
    tfaces = view(patch_to_tfaces,patch)
    connectivity = Vector{Vector{Int32}}(undef, length(tfaces))
    for (k,tface) in enumerate(tfaces)
      face = tface_to_face[tface]
      cell = face_glue.face_to_cell[tface]
      lface = face_glue.face_to_lface[tface]
      pindex = face_glue.cell_to_lface_to_pindex[cell][lface]
      fperm = get_vertex_permutations(Polytope{D-1}(polys[cell],lface))[pindex]

      connectivity[k] = face_to_vertices[face][fperm]
      for i in eachindex(connectivity[k])
        val = get!(glob_to_loc, connectivity[k][i], nv+1)
        connectivity[k][i] = val # map global vertex index to local index
        nv += (val == nv + 1) # increment only if not already present
      end
    end

    vertices = collect(keys(glob_to_loc))
    coords = vertex_coordinates[vertices]

    new_polys[patch], vperm = ReferenceFEs.polytope_from_faces(D,coords,connectivity)
    new_connectivity[patch] = vertices[vperm]
  end

  return new_polys, Table(new_connectivity)
end

function compose_glues(
  glue12::AdaptivityGlue{RefinementGlue,Dc},
  glue23::AdaptivityGlue{RefinementGlue,Dc}
) where Dc
  
  # n2o_faces_map: for each dimension, maps mesh1 faces to mesh3 faces (0 if no match)
  n2o_faces_map = Vector{Vector{Int32}}(undef, Dc+1)
  for d in 0:Dc
    f1_to_f2 = glue12.n2o_faces_map[d+1]
    f2_to_f3 = glue23.n2o_faces_map[d+1]
    f1_to_f3 = zeros(Int32, length(f1_to_f2))
    for (f1, f2) in enumerate(f1_to_f2)
      if !iszero(f2)
        f1_to_f3[f1] = f2_to_f3[f2]
      end
    end
    n2o_faces_map[d+1] = f1_to_f3
  end

  # Child ids: local index of each mesh1 cell within its mesh3 parent
  n3 = length(glue23.o2n_faces_map)
  n2o_cells = n2o_faces_map[Dc+1]
  o2n_faces_map = Arrays.inverse_table(n2o_cells, n3)
  n2o_cell_to_child_id = Arrays.find_local_index(n2o_cells, o2n_faces_map)

  # Placeholder refinement rules (to be replaced later)
  refinement_rules = Fill(WhiteRefinementRule(TRI), n3)
  is_refined = select_refined_cells(n2o_faces_map[end])

  return AdaptivityGlue(
    RefinementGlue(), n2o_faces_map, n2o_cell_to_child_id, 
    refinement_rules, is_refined, o2n_faces_map
  )
end
