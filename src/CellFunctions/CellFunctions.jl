module CellFunctions

export CellFieldValues
export CellBasisValues
export CellPoints

export CellFunction
export CellField
export CellBasis
export CellGeomap

export expand
export varinner
export attachgeomap
export compose

export CellBasisFromSingleInterpolation

using Numa.Helpers
using Numa.FieldValues
using Numa.CellValues
using Numa.CellValues: CellArrayFromUnaryOp
using Numa.CellValues: CellArrayFromBroadcastUnaryOp
using Numa.CellValues: CellArrayFromBoradcastBinaryOp
using Numa.Polynomials

import Base: +, -, *, /, ∘

import Numa.FieldValues: inner, outer
import Numa: evaluate, gradient
import Numa.CellValues: inputcellarray, computesize, computevals!

include("AbstractCellFunctions.jl")
include("Operators.jl")
include("CellBasisFromSingleInterpolation.jl")
include("CellBasisWithGeomap.jl")

end # module CellFunctions
