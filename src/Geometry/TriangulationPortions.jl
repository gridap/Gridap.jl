
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

function get_cell_map(trian::TriangulationPortion)
  cell_map = get_cell_map(trian.oldtrian)
  reindex(cell_map,trian.cell_to_oldcell)
end

function get_cell_id(trian::TriangulationPortion)
  reindex(get_cell_id(trian.oldtrian),trian.cell_to_oldcell)
end

function restrict(f::AbstractArray,trian::TriangulationPortion)
  reindex(f,trian)
end

function get_normal_vector(trian::TriangulationPortion)
  reindex(get_normal_vector(trian.oldtrian),trian.cell_to_oldcell)
end
