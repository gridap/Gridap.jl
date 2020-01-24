
"""
"""
struct TriangulationPortion{Dc,Dp,G} <: Triangulation{Dc,Dp}
  oldtrian::G
  cell_to_oldcell::Vector{Int}
  @doc """
  """
  function TriangulationPortion(oldtrian::Triangulation{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
    new{Dc,Dp,typeof(oldtrian)}(oldtrian,cell_to_oldcell)
  end
end

function get_reffes(trian::TriangulationPortion)
  get_reffes(trian.oldtrian)
end

function get_cell_type(trian::TriangulationPortion)
  reindex(get_cell_type(trian.oldtrian),trian.cell_to_oldcell)
end

function get_cell_coordinates(trian::TriangulationPortion)
  reindex(get_cell_coordinates(trian.oldtrian),trian.cell_to_oldcell)
end
