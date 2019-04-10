module Geometry

export Grid, CartesianGrid, FlexibleUnstructuredGrid
export cells, points, celltypes
export writevtk

import Base: size, getindex, IndexStyle

using StaticArrays: SVector, MVector, @SVector
using WriteVTK
using WriteVTK.VTKCellTypes: VTK_VERTEX
using WriteVTK.VTKCellTypes: VTK_LINE
using WriteVTK.VTKCellTypes: VTK_TRIANGLE
using WriteVTK.VTKCellTypes: VTK_QUAD
using WriteVTK.VTKCellTypes: VTK_TETRA
using WriteVTK.VTKCellTypes: VTK_HEXAHEDRON

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.CellValues: IndexCellValue, IndexCellArray
using Numa.CellFunctions
using Numa.Polytopes

include("Interfaces.jl")
include("CartesianGrids.jl")
include("UnstructuredGrids.jl")
include("Vtkio.jl")

end # module Geometry
