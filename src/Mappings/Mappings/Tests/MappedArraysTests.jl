module MappedArraysTests

using Test
using Gridap.Arrays
using Gridap.Mappings
using Gridap.Mappings: evaluate
using FillArrays
# using Gridap.Arrays: ArrayWithCounter, reset_counter!

a = rand(3,2,4)
test_array(a,a)

a = rand(3,2,4)
a = CartesianIndices(a)
test_array(a,a)

a = rand(3,2)
a = CartesianIndices(a)
c = apply_mapping(-,a)
test_array(c,-a)

a = rand(12)
c = apply_mapping(-,a)
test_array(c,-a)

a = rand(12)
b = rand(12)
c = apply_mapping(-,a,b)
test_array(c,a.-b)

c = apply_mapping(Float64,-,a,b)
test_array(c,a.-b)

a = rand(0)
b = rand(0)
c = apply_mapping(-,a,b)
test_array(c,a.-b)

a = fill(rand(2,3),12)
b = rand(12)
c = apply_mapping(BroadcastMapping(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),0)
b = rand(0)
c = apply_mapping(BroadcastMapping(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),12)
b = rand(12)
c = apply_mapping(BroadcastMapping(-),a,b)
d = apply_mapping(BroadcastMapping(+),a,c)
e = apply_mapping(BroadcastMapping(*),d,c)
test_array(e,[((ai.-bi).+ai).*(ai.-bi) for (ai,bi) in zip(a,b)])

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
v = apply_mapping(a,x,y)
r = [(xi+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)
v = apply_mapping(Float64,a,x,y)
test_array(v,r)

a = Fill(BroadcastMapping(+),10)
x = [rand(2,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply_mapping(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

a = Fill(BroadcastMapping(+),10)
x = [rand(mod(i-1,3)+1,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply_mapping(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

using Gridap.NewFields
using Gridap.NewFields: MockField, MockBasis, OtherMockBasis
using Gridap.TensorValues

np = 4
p = Point(1,2)
x = fill(p,np)

v = 3.0
d = 2
f = MockField{d}(v)
f = MockField(d,v)

a = Fill(FunctionMapping(MockField),5)
b = Fill(2,5)
c = [6.0,2.0,5.0,7.0,9.0]
v = apply_mapping(a,b,c)

@test isa(v,AbstractArray{<:Mapping})

vv = apply_mapping(FunctionMapping(MockField),b,c)

@test vv == v

xx = fill(x,5)
r = apply_mapping(v,xx)
@test all(r[1] .== c[1])

# Test the intermediate results caching mechanism


# a = ArrayWithCounter(fill(rand(2,3),12))
# b = ArrayWithCounter(rand(12))
# c = apply_mapping(bcast(-),a,b)
# d = apply_mapping(bcast(+),a,c)
# e = apply_mapping(bcast(*),d,c)
# r = [ (ai.-bi).*(ai.+(ai.-bi)) for (ai,bi) in zip(a,b)]
# cache = array_cache(e)
# reset_counter!(a)
# reset_counter!(b)
# for i in 1:length(e)
#   ei = getindex!(cache,e,i)
#   ei = getindex!(cache,e,i)
#   ei = getindex!(cache,e,i)
# end

# @test all(a.counter .== 2)
# @test all(b.counter .== 1)

# l = 10
# ai = 1.0
# bi = 2.0
# a = Fill(ai,l)
# b = Fill(bi,l)
# c = apply_mapping(+,a,b)
# r = map(+,a,b)
# test_array(c,r)
# @test isa(c,Fill)

# l = 10
# ai = [8 0; 0 4]
# a = Fill(ai,l)
# c = apply_mapping(inv,a)
# r = map(inv,a)
# test_array(c,r)
# @test isa(c,Fill)

# l = 10
# ai = [8 0; 0 4]
# a = fill(ai,l)
# c = apply_mapping(inv,a)
# r = map(inv,a)
# test_array(c,r)

# l = 0
# ai = [8 0; 0 4]
# a = fill(ai,l)
# c = apply_mapping(inv,a)
# r = map(inv,a)
# test_array(c,r)

# ai = [8 0; 0 4]
# a = CompressedArray([ai,],Int[])
# c = apply_mapping(inv,a)
# r = map(inv,a)
# test_array(c,r)

# l = 10
# ai = [8, 0]
# a = fill(ai,l)
# f(ai) = ai[2]-ai[1]
# c = apply_mapping(f,a)
# r = map(f,a)
# test_array(c,r)

# g(ai) = ai[2]-ai[1]
# import Gridap.Arrays: kernel_testitem!
# kernel_testitem!(c,::typeof(g),ai) = zero(eltype(ai))
# l = 0
# ai = [8, 0]
# a = fill(ai,l)
# c = apply_mapping(g,a)
# r = map(g,a)
# test_array(c,r)


#f = rand(10)
#a = rand(10)
#@test apply_mapping(f,a) === f
#
#f = fill(rand(4),10)
#a = rand(10)
#@test apply_mapping(f,a) === f
#
#f = fill(rand(4),10)
#g = fill(rand(4),10)
#a = rand(10)
#h = apply_all((f,g),a)
#@test h[1] === f
#@test h[2] === g

#l = 10
#ai = 1.0
#bi = 2.0
#a = fill(ai,l)
#b = fill(bi,l)
#c = apply_mapping(+,a,b)
#d = apply_mapping(c,a,b)
#@test c == d
#@test typeof(d) == typeof(c)

end # module
