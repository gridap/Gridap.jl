module MappedArraysTests

using Test
using Gridap.Arrays
using Gridap.Mappings
using Gridap.Mappings: evaluate
using Gridap.Mappings: test_mapped_array
using FillArrays
# using Gridap.Arrays: ArrayWithCounter, reset_counter!

a = rand(3,2,4)
test_array(a,a)

a = rand(3,2,4)
a = CartesianIndices(a)
test_array(a,a)

a = rand(3,2)
a = CartesianIndices(a)
c = apply(-,a)
test_array(c,-a)

a = rand(12)
c = apply(-,a)
test_array(c,-a)

a = rand(12)
b = rand(12)
c = apply(-,a,b)
test_array(c,a.-b)

c = apply(Float64,-,a,b)
test_array(c,a.-b)

a = rand(0)
b = rand(0)
c = apply(-,a,b)
test_array(c,a.-b)

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
c = array_caches(a,b)
i = 1
v = getitems!(c,(a,b),i)
@test c == (nothing,nothing)
@test v == (a[i],b[i])

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
ai = testitem(a)
@test ai == a[1]
ai, bi = testitems(a,b)
@test ai == a[1]
@test bi == b[1]

a = fill(rand(Int,2,3),0)
b = fill(1,0)
ai = testitem(a)
@test ai == Array{Int,2}(undef,0,0)
ai, bi = testitems(a,b)
@test ai == Array{Int,2}(undef,0,0)
@test bi == zero(Int)

a = fill(+,10)
x = rand(10)
y = rand(10)
v = apply(a,x,y)
r = [(xi+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)
v = apply(Float64,a,x,y)
test_array(v,r)

p = 4
p1 = [1,2]
p2 = [2,1]
p3 = [4,3]
p4 = [6,1]

x = [p1,p2,p3,p4]

aa = Fill(fa,4)
r = apply(aa,x)
@test all([ r[i] ≈ 2*x[i] for i in 1:4])

bb = Fill(fb,4)
r = apply(bb,x)
@test all([ r[i] ≈ sqrt.(x[i]) for i in 1:4])

cm = apply(operation,aa,bb)
r = apply(cm,x)
@test all([ r[i] ≈ 2*(sqrt.(x[i])) for i in 1:4])

# BroadcastMapping

a = fill(rand(2,3),12)
b = rand(12)
c = apply(BroadcastMapping(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),0)
b = rand(0)
c = apply(BroadcastMapping(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),12)
b = rand(12)
c = apply(BroadcastMapping(-),a,b)
d = apply(BroadcastMapping(+),a,c)
e = apply(BroadcastMapping(*),d,c)
test_array(e,[((ai.-bi).+ai).*(ai.-bi) for (ai,bi) in zip(a,b)])

a = Fill(BroadcastMapping(+),10)
x = [rand(2,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

a = Fill(BroadcastMapping(+),10)
x = [rand(mod(i-1,3)+1,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

# using Gridap.NewFields
# using Gridap.NewFields: MockField, MockBasis, OtherMockBasis
# using Gridap.TensorValues

# np = 4
# p = Point(1,2)
# x = fill(p,np)

# v = 3.0
# d = 2
# f = MockField{d}(v)
# f = MockField(d,v)

# # a = Fill(FunctionMapping(MockField),5)
# b = Fill(2,5)
# c = [6.0,2.0,5.0,7.0,9.0]
# v = apply_mapping(mock_field,b,c)

# @test isa(v,AbstractArray{<:Mapping})

# # vv = apply_mapping(FunctionMapping(MockField),b,c)
# # @test vv == v

# xx = fill(x,5)
# r = apply(v,xx)
# @test all(r[1] .== c[1])

end # module
