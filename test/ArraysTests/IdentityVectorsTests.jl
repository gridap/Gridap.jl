module IdentityVectorsTests

using Test
using Gridap.Arrays

l = 10
a = IdentityVector(l)
b = collect(1:l)
test_array(a,b)

c = rand(l)
d = lazy_map(Reindex(c),a)
@test d === c

end # module
