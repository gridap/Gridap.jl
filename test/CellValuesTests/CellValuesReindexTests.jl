module CellValuesReindexTests

using Test
using Gridap
using Gridap.CellValuesGallery
using Gridap.CellMapApply: IndexCellMapFromKernel
using ..CellValuesMocks
using ..MapsMocks

a = [1,2,4,5,2,4,6,7,4,3,5,1,2,3]
ca = CellValueFromArray(a)
b = [10,11,3,1]
cb = CellValueFromArray(b)

cc = reindex(ca,cb)

d = a[b]
test_index_cell_value(cc,d)

b = [10,10,10,10]
cb = TestIterCellValue(10,4)

cc = reindex(ca,cb)
d = a[b]
test_iter_cell_value(cc,d)

l = 10
a = VectorValue(10,10)
p1 = VectorValue(1,1)
p2 = VectorValue(2,2)
p3 = VectorValue(3,3)
p = [p1,p2,p3]
m = MockMap(a)
cp = ConstantCellValue(p,l)
cm1 = ConstantCellValue(m,l)
cm2 = CompressedCellValue([m,],fill(1,l))
cm3 = cm1 + cm2
b = [10,11,3,1]
cb = CellValueFromArray(b)
rcm = reindex(cm3,cb)
rcp = reindex(cp,cb)
@test isa(rcm,IndexCellMapFromKernel)

end # module
