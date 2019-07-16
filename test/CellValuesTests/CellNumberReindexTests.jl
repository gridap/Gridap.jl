module CellNumberReindexTests

using Test
using Gridap
using Gridap.CellValuesGallery
using Gridap.CellMapApply: IndexCellMapFromKernel
using ..CellValuesMocks

a = [1,2,4,5,2,4,6,7,4,3,5,1,2,3]
ca = CellValueFromArray(a)
b = [10,11,3,1]
cb = CellValueFromArray(b)

cc = reindex(ca,cb)

d = a[b]
test_index_cell_number(cc,d)

b = [10,10,10,10]
cb = TestIterCellValue(10,4)

cc = reindex(ca,cb)
d = a[b]
test_iter_cell_value(cc,d)

end # module
