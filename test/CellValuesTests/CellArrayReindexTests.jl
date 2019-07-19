module CellArrayReindexTests

using Test
using Gridap
using Gridap.CellValuesGallery
using ..CellValuesMocks

a = [[1,2,4,5],[2,4,6,7],[4,3,5,1],[2,3]]
ca0 = CellValueFromArray(a)
ca = -ca0
b = [3,1,2]
cb = CellValueFromArray(b)

cc = reindex(ca,cb)

d = -a[b]
test_index_cell_array(cc,d)

z = [2,3,1]
cz = CellValueFromArray(z)

cc2 = reindex(ca,cz)

d2 = -a[z]

# This tests that ci and ci2 are properly cached
for (k,(ci,ci2)) in enumerate(zip(cc,cc2))
  @test ci == d[k]
  @test ci2 == d2[k]
end

b = [3,3,3,3]
cb = TestIterCellValue(3,4)

cc = reindex(ca,cb)
d = -a[b]
test_iter_cell_array(cc,d)


end # module
