module CellArrays

using LinearAlgebra: det
import LinearAlgebra
using Base: @propagate_inbounds
using Base.Cartesian: @nloops, @nexprs, @nref

using Numa.Helpers
using Numa.FieldValues

export CellArray
export CellVector
export IndexableCellArray
export CellArrayFromUnaryOp
export CellArrayFromElemUnaryOp
export ConstantCellArray
export cellsize
export celllength
export cellsum
export cellnewaxis
export binner
export bouter

include("Helpers.jl")
include("CachedArray.jl")
include("AbstractCellArrays.jl")
include("Operators.jl")
include("ConstantCellArrays.jl")

end # module CellArrays
