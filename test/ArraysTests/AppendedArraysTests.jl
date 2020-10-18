module AppendedArraysTests

using Test
using Gridap.Arrays

a = collect(Float64,11:20)
b = collect(Float64,101:120)
c = lazy_append(a,b)
r = vcat(a,b)
test_array(c,r)

@test sum(c) == sum(r)

d = lazy_map(-,c)
r = -c
test_array(d,r)
@test isa(d,AppendedArray)

e = lazy_map(-,c,d)
r = c-d
test_array(e,r)
@test isa(e,AppendedArray)

d = lazy_map(-,Float64,c)
r = -c
test_array(d,r)
@test isa(d,AppendedArray)

e = lazy_map(-,Float64,c,d)
r = c-d
test_array(e,r)
@test isa(e,AppendedArray)

a = collect(11:20)
b = collect(Float64,101:120)
c = lazy_append(a,b)
d = lazy_map(i->Float64(i),c)
r = vcat(Float64.(a),b)
test_array(d,r)

end # module
