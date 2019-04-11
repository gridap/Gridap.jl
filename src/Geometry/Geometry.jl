module Geometry

export Grid, CartesianGrid
export UnstructuredGrid, FlexibleUnstructuredGrid
export cells, points, celltypes

import Base: size, getindex, IndexStyle

using StaticArrays: SVector, MVector, @SVector

using Numa
using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.CellValues: IndexCellValue, IndexCellArray
using Numa.CellFunctions
using Numa.Polytopes

include("Interfaces.jl")
include("CartesianGrids.jl")
include("UnstructuredGrids.jl")

end # module Geometry
