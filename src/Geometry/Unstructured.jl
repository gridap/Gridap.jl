module Unstructured

# Dependencies of this module

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.Geometry
using Numa.Geometry.Cartesian

# Functionality provided by this module

export UnstructuredGrid
export FlexibleUnstructuredGrid
import Numa.Geometry: points, cells, celltypes, cellorders

"""
Struct representing an unstructured grid with efficient memory layout
"""
struct UnstructuredGrid{D,Z,A<:IndexCellValue{NTuple{Z,Int}},B<:IndexCellValue{Int}} <: Grid{D,Z}
  points::Vector{Point{D}}
  cells_data::Vector{Int}
  cells_ptrs::Vector{Int}
  ctypes::A
  corders::B
end

points(self::UnstructuredGrid) = CellValueFromArray(self.points)

cells(self::UnstructuredGrid) = CellVectorFromDataAndPtrs(self.cells_data,self.cells_ptrs)

celltypes(self::UnstructuredGrid) = self.ctypes

cellorders(self::UnstructuredGrid) = self.corders

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

"""
Construct an `UnstructuredGrid` from a `CartesianGrid`
"""
function UnstructuredGrid(grid::CartesianGrid{D}) where D
  ps = _compute_points(grid)
  ts = celltypes(grid)
  os = cellorders(grid)
  data, ptrs = _compute_cells(grid,UnstructuredGrid{D,D})
  UnstructuredGrid(ps,data,ptrs,ts,os)
end

"""
Construct an `FlexibleUnstructuredGrid` from a `CartesianGrid`
"""
function FlexibleUnstructuredGrid(grid::CartesianGrid{D}) where D
  ps = _compute_points(grid)
  ts = celltypes(grid)
  os = cellorders(grid)
  cs = _compute_cells(grid,FlexibleUnstructuredGrid{D,D})
  FlexibleUnstructuredGrid(ps,cs,ts,os)
end

# Low level details

function _compute_points(grid::CartesianGrid{D}) where D
  ps = Array{Point{D},1}(undef,(length(points(grid)),))
  for (i,xi) in enumerate(points(grid))
    ps[i] = xi
  end
  ps
end

function _compute_cells(grid::CartesianGrid{D},::Type{FlexibleUnstructuredGrid{D,D}}) where D
  cs = [ Array{Int,1}(undef,(2^D,)) for i in 1:length(cells(grid)) ]
  for (i,ci) in enumerate(cells(grid))
    cs[i] .= ci
  end
  cs
end

function _compute_cells(grid::CartesianGrid{D},::Type{UnstructuredGrid{D,D}}) where D
  ptrs = fill(2^D,(length(cells(grid))+1,))
  length_to_ptrs!(ptrs)
  data = zeros(Int,ptrs[end]-1)
  _fill_cell_data!(data,cells(grid))
  (data, ptrs)
end

function _fill_cell_data!(data,cells)
  k = 1
  for v in cells
    for vi in v
      @inbounds data[k] = vi
      k +=1
    end
  end
end

end # module Unstructured
