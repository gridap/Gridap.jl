module TablesTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Io
using JSON
using FillArrays

data = Float64[2,3,1,3,6,7,3,2,5,6,3,4]
ptrs = [1,4,4,7,13]
a = Table(data,ptrs)
b = [ data[ptrs[i]:ptrs[i+1]-1] for i in 1:length(ptrs)-1]
test_array(a,b)

c = convert(Table{Float64,Vector{Float64},Vector{Int32}},a)
test_array(c,b)

data = Fill(1.3,12)
ptrs = [1,4,4,7,13]
d = Table(data,ptrs)
e = [ data[ptrs[i]:ptrs[i+1]-1] for i in 1:length(ptrs)-1]
test_array(d,e)

data = Int[2,3,1,3,6,7,3,2,5,6,3,4]
ptrs = [1,4,4,7,13]
a = Table(data,ptrs)

a_view = view(a,1)
@test isa(a_view,AbstractVector{Int})
@test length(a_view) == 3
@test a_view == [2,3,1]

a_view = view(a,1:3)
@test isa(a_view,Table{Int})
@test length(a_view) == 3
@test a_view.ptrs == [1,4,4,7]
@test a_view.data == [2,3,1,3,6,7]

data  = reinterpret(Int,Vector{Float64}(undef,12))
data[1:6] .= a.data[7:12]
data[7:9] .= a.data[4:6]
data[10:12] .= a.data[1:3]

perm = Vector{Int}(undef,12)
perm[1:3]  .= 10:12
perm[4:6]  .= 7:9
perm[7:12] .= 1:6
data = lazy_map(Reindex(data),perm)
b = Table(data,ptrs)
test_array(a,b)

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
_ptrs2 = Arrays.generate_ptrs(vv)

data = [1, 2, 3, 2, 3, 5, 8, 1, 2, 4]
ptrs = [1, 4, 6, 8, 8, 11]

@test data == _data
@test ptrs == _ptrs
@test ptrs == _ptrs2

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

a_to_lb_to_b = [[1,2,3],[2,3],[5,8],[2],[1,2,4]]
a_to_lb_to_b = Table(a_to_lb_to_b)
a_to_lb = fill(1,length(a_to_lb_to_b))
a_to_b = lazy_map(getindex,a_to_lb_to_b,a_to_lb)
r = [ lb_to_b[a_to_lb[a]] for (a,lb_to_b) in enumerate(a_to_lb_to_b) ]
test_array(a_to_b,r)

b_to_la_to_a = Table([[5,1,4,2,3],[1,2,3],[5,4],[2,4,5,3],[5,1,2,4]])
a_to_b = [2,5,4,3,1]
a_to_la = find_local_index(a_to_b,b_to_la_to_a)
r = [1,3,4,2,1]
test_array(a_to_la,r)

b_to_la_to_a = Table([[5,1,4,2,3],[1,3,2],[5,4],[2,4,5,3],[5,1,2,4]])
c_to_la_to_a = Table([[1,3],[3,4,5]])
c_to_b = [2,4]
c_to_la_to_lainb = find_local_index(c_to_la_to_a,c_to_b,b_to_la_to_a)
@test c_to_la_to_lainb == [[1,2],[4,2,3]]

a_to_lb_to_b = Table([[5,1,4,2,3],[1,2,3],[5,4],[2,4,5,3],[5,1,2,4]])
a_to_b = [5,2,4,3,1]
a_to_lb = Arrays.find_local_nbor_index(a_to_b,a_to_lb_to_b)
test_array(a_to_lb,[1,2,2,4,2])

a = [[5,1,4,2,3],[1,2,3],[5,4],[2,4,5,3],[5,1,2,4]]
t = Table(a)
r = 2:4
u = t[r]
b = a[r]
@test isa(u,Table)
test_array(u,b)

r = [2,1,3]
u = t[r]
b = a[r]
@test isa(u,Table)
test_array(u,b)

t = empty_table(3)
r = Vector{Int}[[], [], []]
test_array(t,r)

a_to_bs = Table([[1,2,3],[7,8],[4,5,6]])
b_to_a = flatten_partition(a_to_bs)
r = [1, 1, 1, 3, 3, 3, 2, 2]
test_array(b_to_a,r)

data = Float64[2,3,1,3,6,7,3,2,5,6,3,4]
ptrs = [1,4,4,7,13]
a = Table(data,ptrs)

vv = [[1,2,3],[1,2,3]]
a = Table(vv)
b = Arrays.inverse_table(a)
@test b == [[1,2] for i in 1:3]

ptrs = [1,4,8]
x = Arrays.local_identity_array(ptrs)
@test x == [1,2,3,1,2,3,4]
y = Arrays.block_identity_array(ptrs)
@test y == [1,1,1,2,2,2,2]

vv = [[1,2,3],[2,3],[5,8],Int[],[1,2,4]]
a = Table(vv)
@test datarange(a,1) == 1:3
@test datarange(a,1:3) == 1:7
@test dataview(a,1) == [1,2,3]
@test dataview(a,1:3) == [1,2,3,2,3,5,8]
@test collect(dataiterator(a)) == [
  (1,1,1), (1,2,2), (1,3,3), (2,1,2), (2,2,3), (3,1,5), (3,2,8), (5,1,1), (5,2,2), (5,3,4)
]

dict = to_dict(a)
b = from_dict(Table{Float64,Vector{Float64},Vector{Int32}},dict)
@test a == b

s = to_json(a)
b = from_json(Table{Float64,Vector{Float64},Vector{Int32}},s)
@test a == b

d = mktempdir()
f = joinpath(d,"a.jld2")

to_jld2_file(a,f)
@test a == from_jld2_file(typeof(a),f)

# gather_table_values / gather_table_values!

cell_ids = Table([[1,3,5],[2,4,6],[5,6,7]])
cell_vals = [[10.0,30.0,50.0],[20.0,40.0,60.0],[55.0,66.0,70.0]]
vals = gather_table_values(cell_ids, cell_vals, 7)
@test vals[1] == 10.0
@test vals[2] == 20.0
@test vals[3] == 30.0
@test vals[4] == 40.0
@test vals[5] == 55.0  # last-write-wins from cell 3
@test vals[6] == 66.0
@test vals[7] == 70.0

vals2 = zeros(7)
gather_table_values!(vals2, cell_ids, cell_vals)
@test vals2 == vals

# n inferred from data
vals3 = gather_table_values(cell_ids, cell_vals)
@test vals3 == vals

# scatter_table_values (roundtrip)

global_vals = collect(Float64, 100:100:700)
scattered = scatter_table_values(cell_ids, global_vals)
@test scattered[1] == [100.0, 300.0, 500.0]
@test scattered[2] == [200.0, 400.0, 600.0]
@test scattered[3] == [500.0, 600.0, 700.0]

# gather_posneg_table_values / gather_posneg_table_values!

pn_ids = Table([[1,-1,3],[2,-2,4]])
pn_vals = [[10.0,20.0,30.0],[40.0,50.0,60.0]]
(pos, neg) = gather_posneg_table_values(pn_ids, pn_vals)
@test pos == [10.0, 40.0, 30.0, 60.0]
@test neg == [20.0, 50.0]

pos2 = zeros(4); neg2 = zeros(2)
gather_posneg_table_values!(pos2, neg2, pn_ids, pn_vals)
@test pos2 == pos
@test neg2 == neg

# scatter_posneg_table_values (roundtrip)

scattered_pn = scatter_posneg_table_values(pn_ids, pos, neg)
@test scattered_pn[1] == [10.0, 20.0, 30.0]
@test scattered_pn[2] == [40.0, 50.0, 60.0]

# Edge case: empty table

empty_ids = Table(Int[], Int32[1])
empty_vals = [Int[]]
vals_empty = gather_table_values(empty_ids, empty_vals)
@test isempty(vals_empty)

end # module
