module CellValues

export CellValue
export CellArray
export CellVector

export ConstantCellValue
export ConstantCellArray
export ConstantCellVector

export cellsize
export celllength
export cellsum
export cellnewaxis
export binner
export bouter

using Base: @propagate_inbounds
using Base.Cartesian: @nloops, @nexprs, @nref

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

include("Helpers.jl")
include("CachedArray.jl")
include("Interfaces.jl")
include("Operations.jl")
include("ConstantCellValues.jl")

end # module CellValues
