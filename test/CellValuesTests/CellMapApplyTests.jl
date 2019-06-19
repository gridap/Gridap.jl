module CellMapApplyTests

using Test
using Gridap
using Gridap.CachedArrays
using TensorValues

using ..CellValuesMocks
using ..MapsMocks

l = 10

a = VectorValue(10,10)
b = VectorValue(15,20)
p1 = VectorValue(1,1)
p2 = VectorValue(2,2)
p3 = VectorValue(3,3)
p = [p1,p2,p3]
m = MockMap(a)
r = evaluate(m,p)

cm = TestIterCellValue(m,l)
cp = TestIterCellValue(p,l)
rm = [ CachedArray(r) for i in 1:l]

test_iter_cell_map(cm,cp,rm)

cm2 = apply(-,cm,broadcast=true)
rm = [ CachedArray(-r) for i in 1:l]

test_iter_cell_map(cm2,cp,rm)

cm = TestIndexCellValue(m,l)
cp = TestIndexCellValue(p,l)
rm = [ CachedArray(r) for i in 1:l]

test_index_cell_map_with_index_arg(cm,cp,rm)

ca2 = evaluate(cm,cp)

test_index_cell_array(ca2,rm)

cm = TestIterCellValue(m,l)
cp = TestIterCellValue(p,l)
rm = [ CachedArray(r.-p) for i in 1:l]
cm2 = apply(-,cm,cp,broadcast=true)
test_iter_cell_map(cm2,cp,rm)

cm = TestIndexCellValue(m,l)
cp = TestIndexCellValue(p,l)
rm = [ CachedArray(r.-p) for i in 1:l]
cm2 = apply(-,cm,cp,broadcast=true)
test_index_cell_map_with_index_arg(cm2,cp,rm)

a = VectorValue(10.1,10.3)
m = MockMap(a)
r = evaluate(m,p)

cm = TestIterCellValue(m,l)
cp = TestIterCellValue(p,l)
rm = [ CachedArray(r) for i in 1:l]
test_iter_cell_map(cm,cp,rm)

cm2 = apply(-,cm,broadcast=true)
rm = [ CachedArray(-r) for i in 1:l]
test_iter_cell_map(cm2,cp,rm)

cm = TestIndexCellValue(m,l)
cp = TestIndexCellValue(p,l)
rm = [ CachedArray(r) for i in 1:l]
test_index_cell_map(cm,cp,rm)

cm2 = apply(-,cm,broadcast=true)
rm = [ CachedArray(-r) for i in 1:l]
test_index_cell_map_with_index_arg(cm2,cp,rm)

end # module
