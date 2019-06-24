module GridPortions

using Gridap
using Gridap.CellValuesGallery

export GridPortion

import Gridap: points
import Gridap: cells
import Gridap: celltypes
import Gridap: cellorders

struct GridPortion{D,Z} <: Grid{D,Z}
  oldgrid::Grid{D,Z}
  newcell_to_oldcell::IndexCellValue{Int}
  _newcells
  _newcelltypes
  _newcellorders
end

function GridPortion(oldgrid::Grid, newcell_to_oldcell)
  _newcell_to_oldcell = _setup_newcell_to_oldcell(newcell_to_oldcell)
  oldcells = cells(oldgrid)
  oldcelltypes = celltypes(oldgrid)
  oldcellorders = cellorders(oldgrid)
  _newcells = reindex(oldcells, _newcell_to_oldcell)
  _newcelltypes = reindex(oldcelltypes, _newcell_to_oldcell)
  _newcellorders = reindex(oldcellorders, _newcell_to_oldcell)
  GridPortion(
    oldgrid,_newcell_to_oldcell,_newcells,_newcelltypes,_newcellorders)
end

points(g::GridPortion) = points(g.oldgrid)

cells(g::GridPortion) = g._newcells

celltypes(g::GridPortion) = g._newcelltypes

cellorders(g::GridPortion) = g._newcellorders

_setup_newcell_to_oldcell(a::IndexCellValue) = a

_setup_newcell_to_oldcell(a::AbstractVector) = CellValueFromArray(a)

end # module
