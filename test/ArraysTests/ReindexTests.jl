module ReindexTests

using Test
using FillArrays
using Gridap.Arrays
using Gridap.TensorValues

a = [1,2,3]
b = lazy_map(Reindex(a),[3,2,1])
for i=1:3
    b[i]=i
end
@test b == [1,2,3]
@test a == [3,2,1]

a = [[1,2,4,5],[2,4,6,7],[4,3,5,1],[2,3]]
b = [3,1,2]

c = lazy_map(Reindex(a),b)

d = a[b]
test_array(c,d)

a = Fill(30.0,10)
b = [3,1,2]

c = lazy_map(Reindex(a),b)
@test isa(c,Fill)

d = a[b]
test_array(c,d)

a = Fill(1.0,5,5)
b = [13 23; 15 25]

c = lazy_map(Reindex(a),b)
@test isa(c,Fill)

d = a[b]
test_array(c,d)

a = CompressedArray([30,40,10,20,30],[1,2,3,5,3,1,4,2])
b = [3,1,2]

c = lazy_map(Reindex(a),b)
@test isa(c,CompressedArray)

d = a[b]
test_array(c,d)

x = [1,2,3,5,3,1,4,2]
a = lazy_map(-,x)
b = [3,1,2]

c = lazy_map(Reindex(a),b)

@test isa(c,LazyArray)

d = [a[bi] for bi in b]
test_array(c,d)

d = a[b]
test_array(c,d)

# x = [[1,2,4,5],[2,4,6,7],[4,3,5,1],[2,3]]
# a = Table(x)
# b = [3,1,2]
# c = lazy_map(Reindex(a),b)
# d = a[b]
# test_array(c,d)

p = VectorValue(2.0,1.0)
data = [2,3,1,3,4,4,3,2,5,4,3,4]
ptrs = [1,4,4,7,13]
gid_to_val = [p, 2*p, 3*p, -p, p]
lid_to_gid = Table(data,ptrs)
ca = lazy_map(Broadcasting(Reindex(gid_to_val)),lid_to_gid)
a = [ gid_to_val[gids] for gids in lid_to_gid ]
test_array( ca, a )

ids = [3,1,2]
r = lazy_map(Reindex(ca),ids)
@test isa(r,LazyArray)
test_array(r,ca[ids])

i_to_v=[[1, 2], [5, 6], [1, 5], [2, 6], [2, 3], [6, 7], [3, 7], [3, 4], [7, 8], [4, 8], [9, 10], [5, 9], [6, 10], [10, 11], [7, 11], [11, 12], [8, 12], [13, 14], [9, 13], [10, 14], [14, 15], [11, 15], [15, 16], [12, 16]]
j_to_i=Int64[]
lid_to_gid=lazy_map(Reindex(i_to_v),j_to_i)
gid_to_val=VectorValue{2,Float64}[(0.0, 0.25), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (0.0, 0.5), (0.25, 0.5), (0.5, 0.5), (0.75, 0.5), (0.0, 0.75), (0.25, 0.75), (0.5, 0.75), (0.75, 0.75), (0.0, 1.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0)]
f2 = lazy_map(Broadcasting(Reindex(gid_to_val)),lid_to_gid)

test_array(f2,Vector{VectorValue{2,Float64}}[])

values = Float64[]
indices = Int[]
c = lazy_map(Reindex(values),indices)
test_array(c,Float64[])


end # module
