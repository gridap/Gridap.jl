
"""
    struct GridPortion{Dc,Dp,G} <: Grid{Dc,Dp}
      parent_grid::G
      cell_to_parent_cell::Vector{Int32}
      node_to_parent_node::Vector{Int32}
    end
"""
struct GridPortion{Dc,Dp,G} <: Grid{Dc,Dp}
  parent_grid::G
  cell_to_parent_cell::Vector{Int32}
  node_to_parent_node::Vector{Int32}
  cell_to_nodes::Table{Int32,Vector{Int32},Vector{Int32}}
  @doc """
      GridPortion(parent_grid::Grid{Dc,Dp},cell_to_parent_cell::Vector{Int32}) where {Dc,Dp}
  """
  function GridPortion(parent_grid::Grid,cell_to_parent_cell::AbstractVector{<:Integer})

    Dc = num_cell_dims(parent_grid)
    Dp = num_point_dims(parent_grid)

    parent_cell_to_parent_nodes = get_cell_nodes(parent_grid)
    nparent_nodes = num_nodes(parent_grid)
    parent_node_to_coords = get_node_coordinates(parent_grid)

    node_to_parent_node, parent_node_to_node = _find_active_nodes(
      parent_cell_to_parent_nodes,cell_to_parent_cell,nparent_nodes)

    cell_to_nodes = _renumber_cell_nodes(
      parent_cell_to_parent_nodes,parent_node_to_node,cell_to_parent_cell)

    new{Dc,Dp,typeof(parent_grid)}(parent_grid,cell_to_parent_cell,node_to_parent_node,cell_to_nodes)
  end
end

function GridPortion(parent_grid::Grid,parent_cell_to_mask::AbstractArray{Bool})
  cell_to_parent_cell = findall(collect1d(parent_cell_to_mask))
  GridPortion(parent_grid,cell_to_parent_cell)
end

function GridPortion(parent_grid::Grid,parent_cell_to_mask::AbstractVector{Bool})
  cell_to_parent_cell = findall(parent_cell_to_mask)
  GridPortion(parent_grid,cell_to_parent_cell)
end

function OrientationStyle(::Type{GridPortion{Dc,Dp,G}}) where {Dc,Dp,G}
  OrientationStyle(G)
end

function RegularityStyle(::Type{GridPortion{Dc,Dp,G}}) where {Dc,Dp,G}
  RegularityStyle(G)
end

function get_node_coordinates(grid::GridPortion)
  parent_node_to_coords = get_node_coordinates(grid.parent_grid)
  lazy_map(Reindex(parent_node_to_coords),grid.node_to_parent_node)
end

function get_cell_nodes(grid::GridPortion)
  grid.cell_to_nodes
end

function get_reffes(grid::GridPortion)
  get_reffes(grid.parent_grid)
end

function get_cell_type(grid::GridPortion)
  lazy_map(Reindex(get_cell_type(grid.parent_grid)),grid.cell_to_parent_cell)
end

# Helpers

function _find_active_nodes(oldcell_to_oldnodes,cell_to_oldcell,noldnodes)
  oldnode_is_active = fill(false,noldnodes)
  cache = array_cache(oldcell_to_oldnodes)
  for oldcell in cell_to_oldcell
    oldnodes = getindex!(cache,oldcell_to_oldnodes,oldcell)
    for oldnode in oldnodes
      oldnode_is_active[oldnode] = true
    end
  end
  node_to_oldnode = findall(oldnode_is_active)
  oldnode_to_node = fill(UNSET,noldnodes)
  oldnode_to_node[node_to_oldnode] = 1:length(node_to_oldnode)
  (node_to_oldnode, oldnode_to_node)
end

function _renumber_cell_nodes(oldcell_to_oldnodes,oldnode_to_node,cell_to_oldcell)
  ncells = length(cell_to_oldcell)
  cell_to_nodes_ptrs = zeros(Int32,ncells+1)
  cache = array_cache(oldcell_to_oldnodes)
  for (cell,oldcell) in enumerate(cell_to_oldcell)
    oldnodes = getindex!(cache,oldcell_to_oldnodes,oldcell)
    cell_to_nodes_ptrs[cell+1] = length(oldnodes)
  end
  length_to_ptrs!(cell_to_nodes_ptrs)
  ndata = cell_to_nodes_ptrs[end]-1
  cell_to_nodes_data = zeros(Int32,ndata)
  for (cell,oldcell) in enumerate(cell_to_oldcell)
    oldnodes = getindex!(cache,oldcell_to_oldnodes,oldcell)
    a = cell_to_nodes_ptrs[cell]-1
    for (lnode,oldnode) in enumerate(oldnodes)
      node = oldnode_to_node[oldnode]
      @assert node > 0
      cell_to_nodes_data[a+lnode] = node
    end
  end
  Table(cell_to_nodes_data,cell_to_nodes_ptrs)
end

