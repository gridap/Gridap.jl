module Unstructured

# Dependencies of this module

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.Geometry

# Functionality provided by this module

export UnstructuredGrid
export FlexibleUnstructuredGrid
import Numa.Geometry: points, cells, celltypes, cellorders
export cellsdata, cellsptrs

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

cellsdata(self::UnstructuredGrid) = self.cells_data

cellsptrs(self::UnstructuredGrid) = self.cells_ptrs

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


end # module Unstructured
