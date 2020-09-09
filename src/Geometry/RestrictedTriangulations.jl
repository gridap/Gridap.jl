
"""
"""
struct RestrictedTriangulation{Dc,Dp,G} <: Triangulation{Dc,Dp}
  oldtrian::G
  cell_to_oldcell::Vector{Int}
  oldcell_to_cell::Vector{Int}
  void_to_oldcell::Vector{Int}
  memo::Dict

  function RestrictedTriangulation(
    oldtrian::Triangulation{Dc,Dp},
    cell_to_oldcell::Vector{Int},
    oldcell_to_cell::Vector{Int},
    void_to_oldcell::Vector{Int}) where {Dc,Dp}

    new{Dc,Dp,typeof(oldtrian)}(
      oldtrian,cell_to_oldcell,oldcell_to_cell,void_to_oldcell,Dict())
  end

end

get_memo(a::RestrictedTriangulation) = a.memo

function RestrictedTriangulation(
  oldtrian::Triangulation{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
  n_oldcells = num_cells(oldtrian)
  oldcell_to_cell = fill(Int(UNSET),n_oldcells)
  oldcell_to_cell[cell_to_oldcell] .= 1:length(cell_to_oldcell)
  void_to_oldcell = findall(oldcell_to_cell .== UNSET)
  oldcell_to_cell[void_to_oldcell] .= -(1:length(void_to_oldcell))
  RestrictedTriangulation(oldtrian,cell_to_oldcell,oldcell_to_cell,void_to_oldcell)
end

function RestrictedTriangulation(
  oldtrian::Triangulation{Dc,Dp},oldcell_to_mask::Vector{Bool}) where {Dc,Dp}
  cell_to_oldcell = findall(oldcell_to_mask)
  RestrictedTriangulation(oldtrian,cell_to_oldcell)
end

function get_reffes(trian::RestrictedTriangulation)
  get_reffes(trian.oldtrian)
end

function get_cell_type(trian::RestrictedTriangulation)
  reindex(get_cell_type(trian.oldtrian),trian.cell_to_oldcell)
end

function get_cell_coordinates(trian::RestrictedTriangulation)
  reindex(get_cell_coordinates(trian.oldtrian),trian.cell_to_oldcell)
end

function get_cell_id(trian::RestrictedTriangulation)
  trian.cell_to_oldcell
end

function compute_cell_map(trian::RestrictedTriangulation)
  cell_map = get_cell_map(trian.oldtrian)
  ReindexedCellMap(cell_map,trian.cell_to_oldcell)
end

