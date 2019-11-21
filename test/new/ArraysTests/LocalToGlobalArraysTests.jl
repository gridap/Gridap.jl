module LocalToGlobalArraysTests

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

end # module
