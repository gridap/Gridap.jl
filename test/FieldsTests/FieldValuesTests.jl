module FieldValuesTests

using Test
using Gridap
using TensorValues
using LinearAlgebra: dot
using LinearAlgebra: norm

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

t = TensorValue(1,2,3,4)
@test trace(t) == 5
@test tr(t) == 5

t = TensorValue(1,2,3,4,5,6,7,8,9)
@test trace(t) == 15
@test tr(t) == 15

@test symmetic_part(t) == TensorValue(1.0, 3.0, 5.0, 3.0, 5.0, 7.0, 5.0, 7.0, 9.0)

a = TensorValue(1,2,3,4)
b = a'
@test b == TensorValue(1,3,2,4)
@test a*b == TensorValue(10,14,14,20)

u = VectorValue(1.0,2.0)
v = VectorValue(2.0,3.0)
@test dot(u,v) ≈ inner(u,v)
@test norm(u) ≈ sqrt(inner(u,u))

end # module
