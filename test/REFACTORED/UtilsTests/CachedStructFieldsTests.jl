module CachedStructFieldsTest

using Test
using Gridap.CachedStructFields

a = 10

cache = CachedStructField(a)

@test cache.value == 10

cache.value = 20

@test cache.value == 20

end # module
