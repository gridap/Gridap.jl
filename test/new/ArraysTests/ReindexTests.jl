module ReindexTests

using Test
using Gridap.Arrays

a = [[1,2,4,5],[2,4,6,7],[4,3,5,1],[2,3]]
b = [3,1,2]

c = reindex(a,b)

d = a[b]
test_array(c,d)

end # module
