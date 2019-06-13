module CellFieldsOperationsTests

using Gridap
using ..CellFieldsMocks
using ..CellValuesMocks

l = 10
cf = IterCellFieldMock(2,Int,l)

p1 = Point(1,2)
p2 = Point(0,3)
p3 = Point(7,4)
p4 = Point(8,1)
p = [p1,p2,p3,p4]
cp = TestIterCellValue(p,l)

f,_ = iterate(cf)
fg = gradient(f)
r = evaluate(f,p)
rg = evaluate(fg,p)

cf2 = -cf
v = [-r for i in 1:l]
g = [-rg for i in 1:l]
test_iter_cell_field(cf2,cp,v,g)

cf = IndexCellFieldMock(2,Int,l)
cp = TestIndexCellValue(p,l)
cf2 = -cf
#TODO index
test_iter_cell_field(cf2,cp,v,g)

cf = IterCellFieldMock(2,Int,l)
cf2 = cf+cf
v = [r+r for i in 1:l]
g = [rg+rg for i in 1:l]
test_iter_cell_field(cf2,cp,v,g)

cf = IndexCellFieldMock(2,Int,l)
cp = TestIndexCellValue(p,l)
cf2 = cf+cf
#TODO index
test_iter_cell_field(cf2,cp,v,g)

cf = IterCellBasisMock(2,Int,l)
cp = TestIterCellValue(p,l)

f,_ = iterate(cf)
fg = gradient(f)
r = evaluate(f,p)
rg = evaluate(fg,p)
v = [r+r for i in 1:l]
g = [rg+rg for i in 1:l]
cf2 = cf+cf
test_iter_cell_basis(cf2,cp,v,g)

end # module

