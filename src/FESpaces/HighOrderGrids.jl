
"""
    high_order_grid(
      grid::UnstructuredGrid, topo::GridTopology, orders;
      cell_maps=get_cell_map(grid), force_lagrangian_cell_map=true
    )

Create a new `UnstructuredGrid` whose nodes approximate the geometry at a higher
polynomial order. The node positions are obtained by evaluating `cell_maps`
at the reference-space nodes of the high-order Lagrangian reference FEs.

`cell_maps` should be a cell-wise array of fields mapping reference space to
physical space, it defaults to `get_cell_map(grid)`. If `grid` was created
without given `cell_map`, it is probably first order, so suppling a custom
map is needed to produce truly curved high-order geometry.

By default, the resulting geometry is C0 and has polynomial geometric map
of the given `orders`, since a new `cell_map` is computed from the Lagrangian
reffes that define the nodes. If one wishes to keep the `cell_map` given in
argument, set `force_lagrangian_cell_map=false`.

Conforming node numbering is computed via the `ConformingFESpaces.jl` machinery
so that nodes on shared faces (vertices, edges, faces) are not duplicated.
"""
function high_order_grid(
  grid::UnstructuredGrid{Dc,Dp},
  topo::GridTopology,
  orders;
  cell_maps=get_cell_map(grid),
  force_lagrangian_cell_map=true
) where {Dc,Dp}

  @check OrientationStyle(grid) == OrientationStyle(topo)

  old_reffes = get_reffes(grid)
  cell_types  = get_cell_type(grid)
  polytopes   = map(get_polytope, old_reffes)
  new_reffes  = map(p -> LagrangianRefFE(Float64, p, orders), polytopes)

  # Conforming cell-node IDs
  cell_reffe      = expand_cell_data(new_reffes, cell_types)
  cell_conformity = CellConformity(cell_reffe, GradConformity())
  cell_node_ids, num_nodes, _, _, _ = compute_conforming_cell_dofs(
    cell_conformity, topo, FaceLabeling(topo), Int[]
  )

  # Compute node coordinates
  ctype_ref_nodes = map(get_node_coordinates, new_reffes)
  cell_ref_nodes  = expand_cell_data(ctype_ref_nodes, cell_types)
  cell_phys_nodes = lazy_map(evaluate, cell_maps, cell_ref_nodes)
  node_coordinates = gather_table_values(cell_node_ids, cell_phys_nodes, num_nodes)

  if force_lagrangian_cell_map
    UnstructuredGrid(
      node_coordinates, cell_node_ids, new_reffes, cell_types, OrientationStyle(grid)
    )
  else
    UnstructuredGrid(
      node_coordinates, cell_node_ids, new_reffes, cell_types, OrientationStyle(grid),
      cell_maps
    )
  end
end
