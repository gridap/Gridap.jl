module Simplexify

using Gridap
using Gridap.Helpers
using UnstructuredGrids.Kernels: refine_grid_connectivity

export simplexify

function simplexify(grid::Grid{D,Z}) where {D,Z}
  ugrid = UnstructuredGrid(grid)

  ltcell_to_lpoints = _generate_ltcell_to_lpoints(celltypes(ugrid))

  tcells_data, tcells_ptrs = refine_grid_connectivity(
    ugrid.cells_data, ugrid.cells_ptrs, ltcell_to_lpoints)

  order = 1
  _check_order(order,cellorders(ugrid))

  ntcells = length(tcells_ptrs) - 1
  ex = tuple(fill(TET_AXIS,Z)...)
  _ct = ConstantCellValue(ex,ntcells)
  _co = ConstantCellValue(order,ntcells)

  UnstructuredGrid(
    points(ugrid),
    tcells_data,
    tcells_ptrs,
    _ct,
    _co)

end

# Helpers

function _check_order(order,co)
  @notimplemented
end

function _check_order(order,co::ConstantCellValue)
  @notimplementedif co.value != order
end

function _generate_ltcell_to_lpoints(ct)
  @notimplemented
end

function _generate_ltcell_to_lpoints(ct::ConstantCellValue)
  extrusion = ct.value
  if extrusion == (HEX_AXIS, HEX_AXIS)
    ltcell_to_lpoints = [[1,2,3],[4,3,2]]
    return ltcell_to_lpoints
  elseif extrusion == (HEX_AXIS, HEX_AXIS, HEX_AXIS)
    ltcell_to_lpoints = [
      [7,3,2,1], [7,5,2,1], [7,4,3,2], [7,4,8,2], [7,6,5,2], [7,6,8,2]]
  else
    @notimplemented
  end
end

end # module
