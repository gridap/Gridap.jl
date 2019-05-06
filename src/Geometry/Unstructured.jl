module Unstructured

# Dependencies of this module

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.Geometry
using Numa.Polytopes
using UnstructuredGrids

# Functionality provided by this module

export UnstructuredGrid
export FlexibleUnstructuredGrid
import Numa.Geometry: points, cells, celltypes, cellorders
export cellsdata, cellsptrs
export UGrid
import UnstructuredGrids: UGrid

"""
Struct representing an unstructured grid with efficient memory layout
"""
struct UnstructuredGrid{
  D,Z,P<:AbstractVector{Point{D}},A<:IndexCellValue{NTuple{Z,Int}},B<:IndexCellValue{Int}} <: Grid{D,Z}
  points::P
  cells_data::Vector{Int}
  cells_ptrs::Vector{Int}
  ctypes::A
  corders::B
end

function UnstructuredGrid(
  points::AbstractArray{Float64,2},
  cells_data,
  cells_ptrs,
  ctypes,
  corders)

  dim, npoints = size(points)
  k = reshape(points,(dim*npoints,))
  _points = reinterpret(Point{dim},k)

  UnstructuredGrid(
    _points,
    cells_data,
    cells_ptrs,
    ctypes,
    corders)
end

points(self::UnstructuredGrid) = CellValueFromArray(self.points)

cells(self::UnstructuredGrid) = CellVectorFromDataAndPtrs(self.cells_data,self.cells_ptrs)

cellsdata(self::UnstructuredGrid) = self.cells_data

cellsptrs(self::UnstructuredGrid) = self.cells_ptrs

celltypes(self::UnstructuredGrid) = self.ctypes

cellorders(self::UnstructuredGrid) = self.corders


"""
Create a UGrid from a UnstructuredGrid
"""
UGrid(grid::UnstructuredGrid) = _unstructured_grid_to_ugrid(grid)

"""
Struct representing an unstructured grid with efficient memory layout
"""
struct FlexibleUnstructuredGrid{D,Z,A<:IndexCellValue{NTuple{Z,Int}},B<:IndexCellValue{Int}} <: Grid{D,Z}
  points::Vector{Point{D}}
  cells::Vector{Vector{Int}}
  ctypes::A
  corders::B
end

points(self::FlexibleUnstructuredGrid) = CellValueFromArray(self.points)

cells(self::FlexibleUnstructuredGrid) = CellArrayFromArrayOfArrays(self.cells)

celltypes(self::FlexibleUnstructuredGrid) = self.ctypes

cellorders(self::FlexibleUnstructuredGrid) = self.corders

# Helpers

function _unstructured_grid_to_ugrid(grid::UnstructuredGrid{D}) where D

  x = points(grid)
  npoins = length(x)
  coords = reshape(reinterpret(Float64,x),(D,npoins))

  ctypes, refcells = _setup_ctypes_and_refcells(
    celltypes(grid), cellorders(grid))

  UGrid(
    cellsdata(grid),
    cellsptrs(grid),
    ctypes,
    refcells,
    coords)

end

function _setup_ctypes_and_refcells(cell_to_code, cell_to_order)
  @notimplemented
end

function _setup_ctypes_and_refcells(
  cell_to_code::ConstantCellValue,
  cell_to_order::ConstantCellValue)
  order = celldata(cell_to_order)
  @notimplementedif order != 1
  code = celldata(cell_to_code)
  polytope = Polytope(code)
  refcell = RefCell(polytope)
  ncells = length(cell_to_code)
  cell_to_ctype = ConstantCellValue(1,ncells)
  ctype_to_refcell = [refcell]
  (cell_to_ctype, ctype_to_refcell)
end


end # module Unstructured
