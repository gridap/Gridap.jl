module InterfaceTests

using Test
using Gridap.Arrays

a = rand(20,12)

test_array(a,a)
test_array(a,a,≈)

@test array_caches() == ()
@test getitems!((),(),1) == ()

end # module
