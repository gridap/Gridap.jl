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

v = MultiValue{Tuple{1,2}}(10,20)
n = normalvec(v)
@test n == VectorValue(20,-10)
@test meas(v) == sqrt(500)

v = MultiValue{Tuple{2,3}}(1,0,0,1,0,0)
n = normalvec(v)
@test n == VectorValue(0,0,1)
@test meas(v) ≈ 1.0

v = MultiValue{Tuple{2,3}}(1,0,0,1,1,0)
n = normalvec(v)
@test n == VectorValue(-1,0,1)
@test meas(v) ≈ sqrt(2)

end # module
