module ReindexTests

using Test
using FillArrays
using Gridap.Arrays

a = [[1,2,4,5],[2,4,6,7],[4,3,5,1],[2,3]]
b = [3,1,2]

c = reindex(a,b)

d = a[b]
test_array(c,d)

a = Fill(30.0,10)
b = [3,1,2]

c = reindex(a,b)
@test isa(c,Fill)

d = a[b]
test_array(c,d)

a = CompressedArray([30,40,10,20,30],[1,2,3,5,3,1,4,2])
b = [3,1,2]

c = reindex(a,b)
@test isa(c,CompressedArray)

d = a[b]
test_array(c,d)

end # module
