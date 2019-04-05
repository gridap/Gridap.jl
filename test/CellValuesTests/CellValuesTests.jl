module CellValuesTests

using Test
using LinearAlgebra: inv, det

using Numa.FieldValues
using Numa.CellValues

include("Mocks.jl")
include("OperationsTests.jl")
include("ConstantCellValuesTests.jl")

end # module CellValuesTests
