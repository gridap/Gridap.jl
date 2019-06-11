module FieldValuesTests

using Test
using Gridap
using TensorValues

a = 1
@test isa(a,FieldValue)

a = 1.0
@test isa(a,FieldValue)

a = 1.0 - 3.0im
@test isa(a,FieldValue)

a = VectorValue(1,2,3)
@test isa(a,FieldValue)

end # module
