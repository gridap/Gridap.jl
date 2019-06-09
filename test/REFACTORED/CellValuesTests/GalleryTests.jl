module GalleryTests

using Test
using Gridap
using TensorValues
using StaticArrays

v = [ VectorValue(i,i) for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_number( cv, v )

v = [ SVector(i,i) for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_number( cv, v )

v = CartesianIndices((3,10))
cv = CellValueFromArray(v)
test_index_cell_value( cv, v )

v = [ [i,i,i] for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_array( cv, v )

v = [ [i,i,i] for i in 1:10 ]
cv = CellValueFromArray(v)
test_index_cell_array( cv, v )

end # module
