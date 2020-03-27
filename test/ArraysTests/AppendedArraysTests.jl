module AppendedArraysTests

using Test
using Gridap.Arrays

a = collect(11:20)
b = collect(Float64,101:120)
c = lazy_append(a,b)
r = vcat(a,b)
test_array(c,r)

a = collect(Float64,11:20)
b = collect(Float64,101:120)
c = lazy_append(a,b)
r = vcat(a,b)
test_array(c,r)

d = apply(-,c)
r = -c
test_array(d,r)
@test isa(d,AppendedArray)

e = apply(-,c,d)
r = c-d
test_array(e,r)
@test isa(e,AppendedArray)

d = apply(Float64,-,c)
r = -c
test_array(d,r)
@test isa(d,AppendedArray)

e = apply(Float64,-,c,d)
r = c-d
test_array(e,r)
@test isa(e,AppendedArray)

end # module
