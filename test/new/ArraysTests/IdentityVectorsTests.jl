module IdentityVectorsTests

using Test
using Gridap.Arrays

l = 10
a = identity_vector(l)
b = collect(1:l)
test_array(a,b)

c = rand(l)
d = reindex(c,a)
@test d === c

end # module
