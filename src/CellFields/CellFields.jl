module CellFields

export CellFieldValues
export CellBasisValues
export CellPoints
export expand
export inner

using Numa.FieldValues
using Numa.CellArrays

import Numa.FieldValues: inner, outer

include("CellFieldValues.jl")

end # module CellFields
