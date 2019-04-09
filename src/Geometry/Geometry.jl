module Geometry

export Grid, CartesianGrid
export connectivity, celltypes

import Base: size, getindex, IndexStyle

using StaticArrays: SVector, MVector, @SVector

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.CellValues: IndexCellValue, IndexCellArray
using Numa.Polytopes

import Numa: coordinates

include("Interfaces.jl")
include("Cartesian.jl")

end # module Geometry
