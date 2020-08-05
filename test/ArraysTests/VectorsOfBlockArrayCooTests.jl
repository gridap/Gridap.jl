module VectorsOfBlockArrayCooTests

using Test
using Gridap.Arrays
using BlockArrays
using FillArrays
using LinearAlgebra

l = 10
b11 = [  i*[1 2; 3 4]  for i in 1:l ]
b21 = [  i*[5 6; 7 8; 9 10] for i in 1:l]
b32 = [  i*[5 6 7 8; 7 8 1 3] for i in 1:l]
blocks = (b11,b21,b32)
blockids = [(1,1),(2,1),(3,2)]
ax = Fill((blockedrange([2,3,2]), blockedrange([2,4])),l)
al = VectorOfBlockArrayCoo(blocks,blockids,ax)

@test al[Block(1,1)] === b11
@test al[Block(2,1)] === b21
@test al[Block(2),Block(1)] === b21

@test is_zero_block(al,2,2) == true
@test is_zero_block(al,2,1) == false
@test is_nonzero_block(al,2,1) == true
@test is_zero_block(al,Block(1,2)) == true
@test is_zero_block(al,Block(1),Block(2)) == true

@test blocksize(al) == (3,2)
@test blocksize(al,1) == 3
@test blocksize(al,2) == 2

for (i,b) in enumerateblocks(al)
  @test al[i] === b
end

b2 = [  i*[3,4,5,6] for i in 1:l]
blocks = (b2,)
blockids = [(2,)]
ax = Fill((blockedrange([2,4]),),l)
bl = VectorOfBlockArrayCoo(blocks,blockids,ax)

c1 = [  i*[7,8] for i in 1:l]
c2 = [  i*[3,4,5] for i in 1:l]
blocks = (c1,c2)
blockids = [(1,),(2,)]
ax = Fill((blockedrange([2,3,2]),),l)
cl = VectorOfBlockArrayCoo(blocks,blockids,ax)

dl = apply(MulKernel(),al,bl)
test_array(dl,[ a*b for (a,b) in zip(al,bl) ])

dl = apply(MulAddKernel(2,3),al,bl,cl)
test_array(dl,[ 2*a*b + 3*c for (a,b,c) in zip(al,bl,cl) ])

dl = apply(transpose,al)
test_array(dl,transpose.(al))

dl = apply(MulKernel(),al,apply(transpose,al))
test_array(dl,[ a*transpose(a) for a in al ])








end # module
