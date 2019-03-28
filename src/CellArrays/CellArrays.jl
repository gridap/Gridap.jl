module CellArrays

using LinearAlgebra: det
import LinearAlgebra

using Numa.Helpers

export CellArray
export IndexableCellArray
export CellArrayFromUnaryOp
export CellArrayFromElemUnaryOp
export ConstantCellArray
export cellsize
export celllength

include("AbstractCellArrays.jl")

include("Operators.jl")

include("ConstantCellArrays.jl")

end # module CellArrays
