module CachedValuesTest

using Test
using Gridap.CachedValues

a = 10

cache = CachedValue(a)

@test cache.value == 10

cache.value = 20

@test cache.value == 20

end # module
