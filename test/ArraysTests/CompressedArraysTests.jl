module CompressedArraysTests

using Test
using Gridap.Arrays
using FillArrays

values = [10,20,31]
ptrs = [1,2,3,3,2,2]
a = CompressedArray(values,ptrs)
r = values[ptrs]
test_array(a,r)

b = lazy_map(-,a)
test_array(b,-r)
@test isa(b,CompressedArray)

c = lazy_map(*,a,b)
test_array(c,r.*(-r))
@test isa(c,CompressedArray)

c = lazy_map(*,Int,a,b)
test_array(c,r.*(-r))
@test isa(c,CompressedArray)

b = Fill(4,length(ptrs))
c = lazy_map(*,a,b)
test_array(c,a.*b)
@test isa(c,CompressedArray)
@test c.ptrs === a.ptrs

c = lazy_map(*,b,a)
test_array(c,a.*b)
@test isa(c,CompressedArray)
@test c.ptrs === a.ptrs

k = CompressedArray([zero,+,-],copy(ptrs))
r = collect(CompressedArray([0,20,-31],ptrs))
c = lazy_map(evaluate,k,a)
test_array(c,r)
@test isa(c,CompressedArray)

k = CompressedArray([zero,+,-],ptrs)
r = collect(CompressedArray([0,20,-31],ptrs))
c = lazy_map(evaluate,k,a)
test_array(c,r)
@test isa(c,CompressedArray)

values = [10,20]
ptrs = [1,2,1,1,2,2]
b = CompressedArray(values,ptrs)

c = lazy_map(*,a,b)
test_array(c,a.*b)
@test ! isa(c,CompressedArray)

c = lazy_map(*,Int,a,b)
test_array(c,a.*b)
@test ! isa(c,CompressedArray)

c = lazy_map(*,Int,b,a)
test_array(c,a.*b)
@test ! isa(c,CompressedArray)

end # module
