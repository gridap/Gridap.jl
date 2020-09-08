module MappingArraysTests

using Test
using Gridap.Mappings
using FillArrays
using Gridap.TensorValues

using Gridap.NewFields
using Gridap.NewFields: MockField, MockBasis, OtherMockBasis

using Gridap.Arrays

p = 4
p1 = Point(1,2)
p2 = Point(2,1)
p3 = Point(4,3)
p4 = Point(6,1)
x = [p1,p2,p3,p4]

a = FunctionMapping(x -> x[1]^2)
b = FunctionMapping(x -> sqrt(x[1]))
h = composition(a,b)
evaluate(h,x[1])

aa = Fill(a,4)
bb = Fill(b,4)
ff = Fill(FunctionMapping(composition),4)

apply_mapping(ff,aa,bb)
cm = apply_function(composition,aa,bb)
r = apply_mapping(cm,x)

@test all([ r[i] .≈ x[i][1] for i in 1:4])

end #module
