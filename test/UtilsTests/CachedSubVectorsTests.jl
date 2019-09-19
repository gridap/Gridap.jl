module CachedSubVectorsTests

using Test
using Gridap
using Gridap.CachedSubVectors

v = collect(1:10)
c = CachedSubVector(v)
@test v == c

c = CachedSubVector(v,1,10)
@test v == c

locate!(c,3,7)
@test v[3:7] == c

c[1] = -1
@test c[1] == -1

end # module
