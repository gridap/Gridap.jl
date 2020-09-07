
function reindex(trian::Triangulation,indices)
  TriangulationPortion(trian,collect(Int,indices))
end

function reindex(trian::Triangulation,indices::Vector{Int})
  TriangulationPortion(trian,indices)
end

function reindex(trian::Triangulation,indices::IdentityVector)
  trian
end

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
  ReindexedCellMap(cell_map,trian.cell_to_oldcell)
end

function get_cell_id(trian::TriangulationPortion)
  reindex(get_cell_id(trian.oldtrian),trian.cell_to_oldcell)
end

function CellField(object,trian::TriangulationPortion)
  CellField(object,trian.oldtrian)
end

function get_normal_vector(trian::TriangulationPortion)
  nold = get_normal_vector(trian.oldtrian)
  ϕold = get_cell_map(trian.oldtrian)
  ϕ = get_cell_map(trian)
  (nold∘ϕold)∘inverse_map(ϕ)
end
