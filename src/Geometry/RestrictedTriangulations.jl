
"""
"""
struct RestrictedTriangulation{Dc,Dp,G} <: Triangulation{Dc,Dp}
  oldtrian::G
  cell_to_oldcell::Vector{Int}
  oldcell_to_cell::Vector{Int}

  function RestrictedTriangulation(
    oldtrian::Triangulation{Dc,Dp},
    cell_to_oldcell::Vector{Int},
    oldcell_to_cell::Vector{Int}) where {Dc,Dp}
    new{Dc,Dp,typeof(oldtrian)}(oldtrian,cell_to_oldcell,oldcell_to_cell)
  end

end

function RestrictedTriangulation(
  oldtrian::Triangulation{Dc,Dp},cell_to_oldcell::Vector{Int}) where {Dc,Dp}
  n_oldcells = num_cells(oldtrian)
  oldcell_to_cell = fill(Int(UNSET),n_oldcells)
  oldcell_to_cell[cell_to_oldcell] .= 1:length(cell_to_oldcell)
  RestrictedTriangulation(oldtrian,cell_to_oldcell,oldcell_to_cell)
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

function restrict(f::AbstractArray,trian::RestrictedTriangulation)
  reindex(f,trian)
end

function get_cell_id(trian::RestrictedTriangulation)
  trian.cell_to_oldcell
end
