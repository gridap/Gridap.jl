module CellFields

export CellFieldValues
export CellBasisValues
export CellPoints
export expand
export inner

using Numa.Helpers
using Numa.FieldValues
using Numa.CellArrays

import Numa.FieldValues: inner, outer

include("AbstractCellFunctions.jl")
include("Operators.jl")

end # module CellFields
