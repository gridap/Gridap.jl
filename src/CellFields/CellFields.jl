module CellFields

export CellFieldValues
export CellBasisValues
export CellPoints
export expand
export inner

using Numa.Helpers
using Numa.FieldValues
using Numa.CellArrays
using Numa.Polynomials

import Numa.FieldValues: inner, outer

include("AbstractCellFunctions.jl")
include("Operators.jl")
include("CellBasisFromSingleInterpolation.jl")

end # module CellFields
