
"""
    struct GridPortion{Dc,Dp,G} <: Grid{Dc,Dp}
      oldgrid::G
      cell_to_oldcell::Vector{Int}
    end

Renumeration for cells but not for nodes
"""
struct GridPortion{Dc,Dp,G} <: Grid{Dc,Dp}
  oldgrid::G
  cell_to_oldcell::Vector{Int}
  @doc """
  """
  function GridPortion(oldgrid::Grid{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
    new{Dc,Dp,typeof(oldgrid)}(oldgrid,cell_to_oldcell)
  end
end

function OrientationStyle(::Type{GridPortion{Dc,Dp,G}}) where {Dc,Dp,G}
  OrientationStyle(G)
end

function RegularityStyle(::Type{GridPortion{Dc,Dp,G}}) where {Dc,Dp,G}
  RegularityStyle(G)
end

function get_node_coordinates(grid::GridPortion)
  get_node_coordinates(grid.oldgrid)
end

function get_cell_nodes(grid::GridPortion)
  reindex(get_cell_nodes(grid.oldgrid),grid.cell_to_oldcell)
end

function get_reffes(grid::GridPortion)
  get_reffes(grid.oldgrid)
end

function get_cell_type(grid::GridPortion)
  reindex(get_cell_type(grid.oldgrid),grid.cell_to_oldcell)
end

