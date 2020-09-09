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

# a = FunctionMapping(x -> x[1]^2)
# b = FunctionMapping(x -> sqrt(x[1]))
# h = composition(a,b)
# evaluate(h,x[1])

fa(x) = x[1]^2
fb(x) = sqrt(x[1])
fh = composition(fa,fb)
@test evaluate(fh,x[1]) ≈ x[1][1]

aa = Fill(fa,4)
bb = Fill(fb,4)
cm = apply_mapping(composition,aa,bb)
r = apply_mapping(cm,x)
@test all([ r[i] .≈ x[i][1] for i in 1:4])

cm = apply_mapping(composition,aa,bb)
r = apply_mapping(cm,x)
@test all([ r[i] .≈ x[i][1] for i in 1:4])


end #module
