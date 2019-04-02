module CellFunctions

export CellFieldValues
export CellBasisValues
export CellPoints

export CellFunction
export CellField
export CellBasis
export CellGeomap

export expand
export inner
export attachgeomap
export compose

export CellBasisFromSingleInterpolation

using Numa.Helpers
using Numa.FieldValues
using Numa.CellArrays
using Numa.Polynomials

import Numa.FieldValues: inner, outer
import Numa.Polynomials: evaluate, gradient
import Numa.CellArrays: inputcellarray, computesize, computevals!

include("AbstractCellFunctions.jl")
include("Operators.jl")
include("CellBasisFromSingleInterpolation.jl")
include("CellBasisWithGeomap.jl")

end # module CellFunctions
