module LocalToGlobalPosNegArraysTests

using Gridap.Arrays

p = 1
data = [2,3,-1,3,4,4,3,2,5,4,-3,4]
ptrs = [1,4,4,7,13]
gid_to_val_pos = [p, 2*p, 3*p, -p, p]
gid_to_val_neg = [3*p, -p, p]
lid_to_gid = Table(data,ptrs)
ca = LocalToGlobalPosNegArray(lid_to_gid,gid_to_val_pos,gid_to_val_neg)
a = Vector{Int64}[[2, 3, 3], [], [3, -1, -1], [3, 2, 1, -1, 1, -1]]
test_array( ca, a )

end # module
