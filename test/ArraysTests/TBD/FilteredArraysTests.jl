module FilteredArrays

using Test
using Gridap.Arrays
using Gridap
using FillArrays

v1 = collect(1:4)
v2 = collect(5:8)
v3 = collect(9:12)
values = [v1,v2,v3]
ptrs = [1,2,3,3,2,2]
a = CompressedArray(values,ptrs)
r = values[ptrs]
test_array(a,r)

filter = [false, true, true, false]
filters = Fill(filter,6)

fa = FilteredCellArray(a,filters)

test_array(fa,[a[i][2:3] for i in 1:6])

k = FilterKernel()

cache = kernel_cache(k,filters[1],a[1])

kernel_return_type(k,filters[1],a[1])

lazy_map_kernel!(cache,k,filters[1],a[1])

test_kernel(k,(filters[1],a[1]),a[1][2:3])

@test lazy_map(k,filters,a) == fa

#

v1 = [ [2  4  3  3]; [2  4  3  3]]
v2 = [ [5  8  4  3]; [2  4  4  4]]
v3 = [ [2  2  8  8]; [2  9  8  6]]

values = [v1,v2,v3]
ptrs = [1,2,3,3,2,2]
a = CompressedArray(values,ptrs)
r = values[ptrs]
test_array(a,r)

b = Array{Bool,2}(undef,2,4)
b .= true
b[2,3] = b[1,4] = false

filters = Fill(b,6)

fa = FilteredCellArray(a,filters)

r1 = [2, 2, 4, 4, 3, 3]
r2 = [5, 2, 8, 4, 4, 4]
r3 = [2, 2, 2, 9, 8, 6]
res = [r1,r2,r3]
r = CompressedArray(res,ptrs)
test_array(fa,r)

@test lazy_map(k,filters,a) == fa

end #module
