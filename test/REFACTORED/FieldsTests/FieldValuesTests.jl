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

p = Point(1,2,3)

@test isa(p,VectorValue{3,Int})

p = Point{3,Float64}(1,2,3)

@test isa(p,VectorValue{3,Float64})


end # module
