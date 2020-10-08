module MappingsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using LinearAlgebra
using BlockArrays

test_mapping(+,(3,2),5)

@test return_types((+,/),1,1) == (Int,Float64)

f = bcast(+)
a = rand(3,2)
b = 3
c = a .+ b
test_mapping(f,(a,b),c)

#test_mapping(1,(),1)
#test_mapping(1,(1,),1)
#
#test_mapping([1,2],(),[1,2])
#test_mapping([1,2],(1,),[1,2])

k = elem(-)
test_mapping(k,(1,),-1)
test_mapping(k,([1,2],),[-1,-2])
test_mapping(k,(1,2),-1)
test_mapping(k,(1.0,2),-1.0)
test_mapping(k,(1,2.0),-1.0)
test_mapping(k,([1,2],2),[-1,0])
test_mapping(k,(2,[1,2]),[1,0])
test_mapping(k,([3,4],[1,2]),[2,2])

k = contract(-)
test_mapping(k,([3,4],[1,2]),3-1+4-2)

f = bcast(⋅)
a = fill(TensorValue(2,0,0,0,2,0,0,0,2),2)
b = VectorValue(1,2,3)
c = zeros(VectorValue{3,Int},2)
broadcast!(⋅,c,a,b)
test_mapping(f,(a,b),c)

a = rand(3,4)
b = rand(4)
c = rand(3)
k = MulMapping()
test_mapping(k,(a,b),a*b)
k = MulAddMapping(2,3)
test_mapping(k,(a,b,c),2*a*b+3*c,≈)

a = rand(3,4)
b = rand(4)
k = MulMapping()
test_mapping(k,(a,b),a*b)

blocks = [ [1 2; 3 4], [5 6 7 8; 9 10 11 12; 13 14 15 16], [1 2 3 4; 5 6 7 8], [1 2 3; 4 5 6; 7 8 9] ]
blockids = [(1,1),(2,2),(1,2),(3,3)]
ax = (blockedrange([2,3,3]), blockedrange([2,4,3]))
a = BlockArrayCoo(blocks,blockids,ax)

blocks = [ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
axs = (blockedrange([2,4,3]),)
b = BlockArrayCoo(blocks,blockids,axs)
test_mapping(k,(a,b),a*b)

c = a*b
k = MulAddMapping(2,3)
test_mapping(k,(a,b,c),2*a*b+3*c)

#k = MulMapping()
#cache = return_cache(k,a,b)
#using BenchmarkTools
#@btime evaluate!($cache,$k,$a,$b)
#
#
#k = MulAddMapping(2,3)
#cache = return_cache(k,a,b,c)
#@btime evaluate!($cache,$k,$a,$b,$c)



end # module
