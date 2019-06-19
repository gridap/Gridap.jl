module ConstantCellFieldsTests

using Gridap

using ..FieldsMocks

l = 10
f = MockField(2,Float64)

cf = ConstantCellValue(f,l)

p1 = Point(1,2)
p2 = Point(0,3)
p3 = Point(7,4)
p4 = Point(8,1)
p = [p1,p2,p3,p4]

cp = ConstantCellValue(p,l)

r = evaluate(f,p)
rg = evaluate(∇(f),p)
v = [r for i in 1:l]
g = [rg for i in 1:l]

test_index_cell_field(cf,cp,v,g)

end # module
