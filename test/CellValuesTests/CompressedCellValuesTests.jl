include("../../src/CellValues/CompressedCellValues.jl")

module CompressedCellValuesTests

using Test
using Gridap
using ..CompressedCellValues

values = [10,20,31]
ptrs = [1,2,3,3,2,2]
r = values[ptrs]
cv = IterCompressedCellValue(values,ptrs)
test_iter_cell_value(cv,r)

end # module
