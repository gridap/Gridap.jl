module TablesTests

using Test
using Gridap.Arrays

data = Float64[2,3,1,3,6,7,3,2,5,6,3,4]
ptrs = [1,4,4,7,13]
a = Table(data,ptrs)
b = [ data[ptrs[i]:ptrs[i+1]-1] for i in 1:length(ptrs)-1]
test_array(a,b)

c = convert(Table{Float64,Int32},a)
test_array(c,b)

vv = Array{Array{Int,2},2}(undef,2,2)
vv[1,1] = rand(1:10,2,2)
vv[2,1] = rand(1:10,2,2)
vv[2,2] = rand(1:10,2,2)
vv[1,2] = rand(1:10,2,2)
b = [[ vv[i][j] for j in 1:length(vv[i]) ] for i in 1:length(vv)]
a = Table(vv)
test_array(a,b)

data = Float64[2,3,1,3,6,7,3,2,5,6,3,4]
ptrs = Int32[1,4,4,7,13]
a = Table(data,ptrs)
b = [ data[ptrs[i]:ptrs[i+1]-1] for i in 1:length(ptrs)-1]
test_array(a,b)

data = Float64[]
ptrs = [1,]
a = Table(data,ptrs)
b = [ data[ptrs[i]:ptrs[i+1]-1] for i in 1:length(ptrs)-1]
test_array(a,b)

vv = [[1,2,3],[2,3],[5,8],Int[],[1,2,4]]
a = Table(vv)
test_array(a,vv)

_data, _ptrs = generate_data_and_ptrs(vv)

data = [1, 2, 3, 2, 3, 5, 8, 1, 2, 4]
ptrs = [1, 4, 6, 8, 8, 11]

@test data == _data
@test ptrs == _ptrs

a = [9,2,1,2,4,7,4]
b = [1,9,2,1,2,4,7]

rewind_ptrs!(a)
@test a == b

a = [3,2,4,2]
b = [1,3,7,9]

length_to_ptrs!(a)
@test a == b

pa = [1,3,5,7,9]
pb = [1,3,5,7]

pc = append_ptrs(pa,pb)

@test pa == [1,3,5,7,9]
@test pb == [1,3,5,7]
@test pc == [1, 3, 5, 7, 9, 11, 13, 15]

vv1 = [[1,2,3],[2,3],[5,8],Int[],[1,2,4]]
table1 = Table(vv1)
vv2 = [[1,3],[4,2,3],Int[],Int[],[1,2,4]]
table2 = Table(vv2)

table3 = append_tables_globally(table1,table2)
@test table3 == vcat(table1,table2)

table4 = append_tables_locally((0,5),(table1,table2))
@test table4 == [[1, 2, 3, 6, 8], [2, 3, 9, 7, 8], [5, 8], Int[], [1, 2, 4, 6, 7, 9]]

end # module
