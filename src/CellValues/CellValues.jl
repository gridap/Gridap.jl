module CellValues

export CellValue
export CellArray
export CellVector

export ConstantCellValue
export ConstantCellArray
export ConstantCellVector

export cellsize
export celllength

using Numa.Helpers

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle
import Base: +, -, *, /
import Base: ==
import LinearAlgebra: inv, det

import Numa.FieldValues: inner, outer

include("Interfaces.jl")
include("ConstantCellValues.jl")
include("Operations.jl")

end # module CellValues
