module OtherCellArrays

using LinearAlgebra: det
import LinearAlgebra

using Numa.Helpers

export OtherCellArray
export IndexableCellArray
export OtherCellArrayFromUnaryOp
export OtherCellArrayFromElemUnaryOp
export OtherConstantCellArray
export maxsize
export maxlength

include("AbstractCellArrays.jl")

include("Operators.jl")

include("ConstantCellArrays.jl")

end # module OtherCellArrays
