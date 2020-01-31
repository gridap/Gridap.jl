
"""
    struct GridPortion{Dc,Dp,G} <: Grid{Dc,Dp}
      oldgrid::G
      cell_to_oldcell::Vector{Int}
      node_to_oldnode::Vector{Int}
    end
"""
struct GridPortion{Dc,Dp,G} <: Grid{Dc,Dp}
  oldgrid::G
  cell_to_oldcell::Vector{Int}
  node_to_oldnode::Vector{Int}
  cell_to_nodes::Table{Int,Int32}
  @doc """
      GridPortion(oldgrid::Grid{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
  """
  function GridPortion(oldgrid::Grid{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}

    oldcell_to_oldnodes = get_cell_nodes(oldgrid)
    noldnodes = num_nodes(oldgrid)
    oldnode_to_coords = get_node_coordinates(oldgrid)

    node_to_oldnode, oldnode_to_node = _find_active_nodes(oldcell_to_oldnodes,cell_to_oldcell,noldnodes)
    cell_to_nodes = _renumber_cell_nodes(oldcell_to_oldnodes,oldnode_to_node,cell_to_oldcell)

    new{Dc,Dp,typeof(oldgrid)}(oldgrid,cell_to_oldcell,node_to_oldnode,cell_to_nodes)
  end
end

function OrientationStyle(::Type{GridPortion{Dc,Dp,G}}) where {Dc,Dp,G}
  OrientationStyle(G)
end

function RegularityStyle(::Type{GridPortion{Dc,Dp,G}}) where {Dc,Dp,G}
  RegularityStyle(G)
end

function get_node_coordinates(grid::GridPortion)
  oldnode_to_coords = get_node_coordinates(grid.oldgrid)
  reindex(oldnode_to_coords,grid.node_to_oldnode)
end

function get_cell_nodes(grid::GridPortion)
  grid.cell_to_nodes
end

function get_reffes(grid::GridPortion)
  get_reffes(grid.oldgrid)
end

function get_cell_type(grid::GridPortion)
  reindex(get_cell_type(grid.oldgrid),grid.cell_to_oldcell)
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
  cell_to_nodes_data = zeros(Int,ndata)
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

