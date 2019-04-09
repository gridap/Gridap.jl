module Geometry

import Base: size, getindex, IndexStyle

using StaticArrays: SVector, MVector, @SVector

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues

include("Interfaces.jl")
include("Cartesian.jl")

end # module Geometry
