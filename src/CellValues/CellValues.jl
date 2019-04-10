module CellValues

export CellValue
export CellArray
export CellMatrix
export CellVector
export CellData
export IterData

export ConstantCellValue
export ConstantCellArray
export ConstantCellVector
export ConstantCellMatrix
export ConstantCellData

export CellValueFromArray
export CellArrayFromArrayOfArrays
export CellVectorFromDataAndPtrs

export apply
export celldata
export cellsize
export celllength
export cellsum
export cellnewaxis
export cellmean

using Base: @propagate_inbounds
using Base.Cartesian: @nloops, @nexprs, @nref

using Numa.Helpers

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex, setindex!
import Base: IndexStyle, IteratorSize
import Base: +, -, *, /
import Base: ==
import LinearAlgebra: inv, det

import Numa: flatten
import Numa.FieldValues: inner, outer, meas

include("Helpers.jl")
include("CachedArray.jl")
include("Interfaces.jl")
include("Operations.jl")
include("ConstantCellValues.jl")
include("Wrappers.jl")

end # module CellValues
