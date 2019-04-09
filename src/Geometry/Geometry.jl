module Geometry

import Base: size, getindex, IndexStyle

using StaticArrays: SVector, MVector, @SVector

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.CellValues: IndexCellValue

include("Interfaces.jl")
include("Cartesian.jl")

end # module Geometry
