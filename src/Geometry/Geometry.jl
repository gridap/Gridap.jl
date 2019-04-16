module Geometry

export Grid, CartesianGrid
export UnstructuredGrid, FlexibleUnstructuredGrid
export cells, points, celltypes, gridgraph
export GridGraph
export celltovefs
export veftocells

import Base: size, getindex, IndexStyle

using StaticArrays: SVector, MVector, @SVector

using Numa
using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.CellValues: IndexCellValue, IndexCellArray, IndexCellVector
using Numa.CellFunctions
using Numa.Polytopes
using Numa.Meshes

include("Interfaces.jl")
include("CartesianGrids.jl")
include("UnstructuredGrids.jl")

end # module Geometry
