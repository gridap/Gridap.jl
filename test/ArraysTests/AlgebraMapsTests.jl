module AlgebraMapsTests

using Test
using Gridap.Arrays

a = rand(3,4)
b = rand(4,5)

c = a*b
test_mapping(c,*,a,b)
cache = return_cache(*,a,b)
@test evaluate!(cache,*,a,b) === evaluate!(cache,*,a,b)

a = rand(3,4)
b = rand(4)

c = a*b
test_mapping(c,*,a,b)
cache = return_cache(*,a,b)
@test evaluate!(cache,*,a,b) === evaluate!(cache,*,a,b)

a = 3
b = 4

c = a*b
test_mapping(c,*,a,b)

end # module
