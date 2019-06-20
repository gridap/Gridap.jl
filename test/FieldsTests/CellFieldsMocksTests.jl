module CellFieldsMocksTests

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
v = [r for i in 1:l]
g = [rg for i in 1:l]

test_iter_cell_field(cf,cp,v,g)
test_iter_cell_field_without_grad(cf,cp,v)

cf = IndexCellFieldMock(2,Int,l)
cp = TestIndexCellValue(p,l)
test_index_cell_field(cf,cp,v,g)
test_index_cell_field_without_grad(cf,cp,v)

cf = IterCellBasisMock(2,Int,l)
cp = TestIterCellValue(p,l)

f,_ = iterate(cf)
fg = gradient(f)
r = evaluate(f,p)
rg = evaluate(fg,p)
v = [r for i in 1:l]
g = [rg for i in 1:l]

test_iter_cell_basis(cf,cp,v,g)
test_iter_cell_basis_without_grad(cf,cp,v)

cf = IndexCellBasisMock(2,Int,l)
cp = TestIndexCellValue(p,l)
test_index_cell_basis(cf,cp,v,g)
test_index_cell_basis_without_grad(cf,cp,v)

cf = IterCellGeomapMock(2,Int,l)
cp = TestIndexCellValue(p,l)

f,_ = iterate(cf)
fg = gradient(f)
r = evaluate(f,p)
rg = evaluate(fg,p)
v = [r for i in 1:l]
g = [rg for i in 1:l]

test_iter_cell_field(cf,cp,v,g)

cf = IndexCellGeomapMock(2,Int,l)
test_index_cell_field(cf,cp,v,g)

end # module
