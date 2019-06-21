module CellMapsTests

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
@test !isa(cm,CellNumber)
@test !isa(cm,IterCellNumber)
@test !isa(cm,CellArray)
@test !isa(cm,IterCellArray)

cm = TestIndexCellValue(m,l)
cp = TestIndexCellValue(p,l)
rm = [ CachedArray(r) for i in 1:l]
test_index_cell_map_with_index_arg(cm,cp,rm)
test_iter_cell_map_with_index_result(cm,cp,rm)

end # module
