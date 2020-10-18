module PosNegReindexTests

using Test
using Gridap.Arrays
using Gridap.TensorValues

indices = [1,3,-2,2,-1]
values_pos = Float64[40,30,10]
values_neg = -Float64[40,30]
c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices ]
test_array(c,r)

subindices = [1,4,2]
d = lazy_map(Reindex(c),subindices)
@test d === values_pos

subindices = [5,3]
d = lazy_map(Reindex(c),subindices)
@test d === values_neg

subindices = [4,2]
d = lazy_map(Reindex(c),subindices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
test_array(d,r)

subindices = [5,]
d = lazy_map(Reindex(c),subindices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
test_array(d,r)

subindices = [1,5,3,4]
d = lazy_map(Reindex(c),subindices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
test_array(d,r)

# Now with PosNegPartition

indices_pos = [1,4,2]
indices = PosNegPartition(indices_pos,5)
test_array(indices,collect(indices))
c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices ]
test_array(c,r)

# This should be for free thanks to the PosNegPartition
d = lazy_map(Reindex(c),indices.ipos_to_i)
@test d === values_pos

# This should be for free thanks to the PosNegPartition
d = lazy_map(Reindex(c),indices.ineg_to_i)
@test d === values_neg

subindices = [4,2]
d = lazy_map(Reindex(c),subindices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
test_array(d,r)

subindices = [5,]
d = lazy_map(Reindex(c),subindices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
test_array(d,r)

subindices = [1,5,3,4]
d = lazy_map(Reindex(c),subindices)
r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
test_array(d,r)

p = 1
data = [2,3,-1,3,4,4,3,2,5,4,-3,4]
ptrs = [1,4,4,7,13]
gid_to_val_pos = [p, 2*p, 3*p, -p, p]
gid_to_val_neg = [3*p, -p, p]
lid_to_gid = Table(data,ptrs)
ca = lazy_map(Broadcasting(PosNegReindex(gid_to_val_pos,gid_to_val_neg)),lid_to_gid)
a = Vector{Int64}[[2, 3, 3], [], [3, -1, -1], [3, 2, 1, -1, 1, -1]]
test_array( ca, a )

i_to_v=[[1, 2], [5, 6], [1, 5], [2, 6], [2, 3], [6, 7], [3, 7], [3, 4], [7, 8], [4, 8], [9, 10], [5, 9], [6, 10], [10, 11], [7, 11], [11, 12], [8, 12], [13, 14], [9, 13], [10, 14], [14, 15], [11, 15], [15, 16], [12, 16]]
j_to_i=Int64[]
lid_to_gid=lazy_map(Reindex(i_to_v),j_to_i)
gid_to_val=VectorValue{2,Float64}[(0.0, 0.25), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (0.0, 0.5), (0.25, 0.5), (0.5, 0.5), (0.75, 0.5), (0.0, 0.75), (0.25, 0.75), (0.5, 0.75), (0.75, 0.75), (0.0, 1.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0)]
f2=lazy_map(Broadcasting(PosNegReindex(gid_to_val,gid_to_val)),lid_to_gid)
test_array(f2,Vector{VectorValue{2,Float64}}[])

indices = Int[]
values_pos = Float64[]
values_neg = Float64[]
c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
test_array(c,Float64[])

## Testing some cases where PosNegReindex can be type-instable
#
#indices = [1,3,-2,2,-1]
#values_pos = Int[40,30,10]
#values_neg = -Float64[40,30]
#c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
#d = lazy_map(i->Float64(i),c)
#display(d)
#
#d = lazy_map(*,c,ones(Float64,size(c)))
#display(d)
#
#indices = [1,1,-1,1,-1]
#values_pos = [+]
#values_neg = [-]
#c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
#d = lazy_map(evaluate, c, rand(size(c)), rand(size(c)))
#display(d)

end # module
