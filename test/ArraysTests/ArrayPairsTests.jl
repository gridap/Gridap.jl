module ArrayPairsTests

using Test
using Gridap.Arrays

a = collect(1:10)
b = collect(11:20)
r = [(i,i+10) for i in 1:10]

c = pair_arrays(a,b)
test_array(c,r)

_a, _b = unpair_arrays(r)
test_array(_a,a)
test_array(_b,b)

_a, _b = unpair_arrays(c)
@test a === _a
@test b === _b

end # module
