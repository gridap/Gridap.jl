
"""
    high_order_grid(
      grid::UnstructuredGrid, topo::GridTopology, orders; 
      cell_maps=get_cell_map(grid)
    )

Create a new `UnstructuredGrid` whose geometry is represented at a higher
polynomial order. The node positions are obtained by evaluating `cell_maps`
(defaulting to the grid's own cell map) at the reference-space nodes of the
high-order Lagrangian reference FEs. 

`cell_maps` should be a cell-wise array of fields mapping reference space to physical
space. Defaults to `get_cell_map(grid)`, which will generally be low-order. 
Supply a custom map to produce truly curved high-order geometry. Note that the arbitrarily 
smooth maps can be provided, but the resulting geometry will only have the polynomial
order of the reference FE nodes. 

Conforming node numbering is computed via the `ConformingFESpaces.jl` machinery
so that nodes on shared faces (vertices, edges, faces) are not duplicated.
"""
function high_order_grid(
  grid::UnstructuredGrid{Dc,Dp},
  topo::GridTopology,
  orders;
  cell_maps=get_cell_map(grid)
) where {Dc,Dp}

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

  UnstructuredGrid(
    node_coordinates, cell_node_ids, new_reffes, cell_types, OrientationStyle(grid)
  )
end
