module LocalToGlobalArraysTests

using Test
using Gridap.TensorValues
using Gridap.Arrays

p = VectorValue(2.0,1.0)
data = [2,3,1,3,4,4,3,2,5,4,3,4]
ptrs = [1,4,4,7,13]
gid_to_val = [p, 2*p, 3*p, -p, p]
lid_to_gid = Table(data,ptrs)
ca = LocalToGlobalArray(lid_to_gid,gid_to_val)
a = [ gid_to_val[gids] for gids in lid_to_gid ]
test_array( ca, a )

ids = [3,1,2]
r = reindex(ca,ids)
@test isa(r,LocalToGlobalArray)
test_array(r,ca[ids])

i_to_v=[[1, 2], [5, 6], [1, 5], [2, 6], [2, 3], [6, 7], [3, 7], [3, 4], [7, 8], [4, 8], [9, 10], [5, 9], [6, 10], [10, 11], [7, 11], [11, 12], [8, 12], [13, 14], [9, 13], [10, 14], [14, 15], [11, 15], [15, 16], [12, 16]]
j_to_i=Int64[]
lid_to_gid=reindex(i_to_v,j_to_i)
gid_to_val=VectorValue{2,Float64}[(0.0, 0.25), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (0.0, 0.5), (0.25, 0.5), (0.5, 0.5), (0.75, 0.5), (0.0, 0.75), (0.25, 0.75), (0.5, 0.75), (0.75, 0.75), (0.0, 1.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0)]
f2=LocalToGlobalArray(lid_to_gid,gid_to_val)
test_array(f2,Vector{VectorValue{2,Float64}}[])

end # module
