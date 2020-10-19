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

#@fverdugo to be deleted since it is kind of type-piracy
#@test al[Block(1,1)] === b11
#@test al[Block(2,1)] === b21
#@test al[Block(2),Block(1)] === b21

#@test is_zero_block(al,2,2) == true
#@test is_zero_block(al,2,1) == false
#@test is_nonzero_block(al,2,1) == true
#@test is_zero_block(al,Block(1,2)) == true
#@test is_zero_block(al,Block(1),Block(2)) == true
#
#@test blocksize(al) == (3,2)
#@test blocksize(al,1) == 3
#@test blocksize(al,2) == 2
#
#for (i,b) in enumerateblocks(al)
#  @test al[i] === b
#end

b1 = [  i*[3,4] for i in 1:l]
b2 = [  i*[3,4,5,6] for i in 1:l]
blocks = (b1,b2,)
blockids = [(1,),(2,)]
ax = Fill((blockedrange([2,4]),),l)
bl = VectorOfBlockArrayCoo(blocks,blockids,ax)

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

dl = lazy_map(*,al,bl)
test_array(dl,[ a*b for (a,b) in zip(al,bl) ])

dl = lazy_map(MulAddMap(2,3),al,bl,cl)
test_array(dl,[ 2*a*b + 3*c for (a,b,c) in zip(al,bl,cl) ])

dl = lazy_map(transpose,al)
test_array(dl,transpose.(al))

dl = lazy_map(*,al,lazy_map(transpose,al))
test_array(dl,[ a*transpose(a) for a in al ])
#print_op_tree(dl)

# in-homogeneous case

l1 = 4
l2 = 6
a11 = vcat([i*[1 2; 3 4] for i in 1:l1 ],[i*[1 2 3; 3 4 5; 6 7 8] for i in 1:l2])
a12 = vcat([10*i*[1 2; 3 4] for i in 1:l1 ],[10*i*[1 2 3; 3 4 5; 6 7 8] for i in 1:l2])
blocks = (a11,a12)
blockids = [(1,1),(1,2)]
ax1 = (blockedrange([2,2]),blockedrange([2,2]))
ax2 = (blockedrange([3,3]),blockedrange([3,3]))
ax = CompressedArray([ax1,ax2],vcat(fill(1,l1),fill(2,l2)))
al = VectorOfBlockArrayCoo(blocks,blockids,ax)

b2 = vcat([i*[1,2] for i in 1:l1 ],[i*[1,2,3] for i in 1:l2])
blocks = (b2,)
blockids = [(2,)]
ax1 = (blockedrange([2,2]),)
ax2 = (blockedrange([3,3]),)
ax = CompressedArray([ax1,ax2],vcat(fill(1,l1),fill(2,l2)))
bl = VectorOfBlockArrayCoo(blocks,blockids,ax)

dl = lazy_map(*,al,bl)
test_array(dl,[ a*b for (a,b) in zip(al,bl) ])

dl = lazy_map(MulAddMap(2,3),al,bl,bl)
test_array(dl,[ 2*a*b + 3*c for (a,b,c) in zip(al,bl,bl) ])

dl = lazy_map(transpose,al)
test_array(dl,transpose.(al))

dl = lazy_map(*,al,lazy_map(transpose,al))
test_array(dl,[ a*transpose(a) for a in al ])

# Blocks of blocks (in-homogeneous case)

blocks = (al,)
blockids = [(1,2),]
_ax1 = (blockedrange([2,2]),blockedrange([2,2]))
_ax2 = (blockedrange([3,3]),blockedrange([3,3]))
ax1 = (blockedrange([_ax1[1],_ax1[1]]),blockedrange([_ax1[2],_ax1[2]]))
ax2 = (blockedrange([_ax2[1],_ax2[1]]),blockedrange([_ax2[2],_ax2[2]]))
ax = CompressedArray([ax1,ax2],vcat(fill(1,l1),fill(2,l2)))
aBl = VectorOfBlockArrayCoo(blocks,blockids,ax)

blocks = (bl,)
blockids = [(1,)]
_ax1 = (blockedrange([2,2]),)
_ax2 = (blockedrange([3,3]),)
ax1 = (blockedrange([_ax1[1],_ax1[1]]),)
ax2 = (blockedrange([_ax2[1],_ax2[1]]),)
ax = CompressedArray([ax1,ax2],vcat(fill(1,l1),fill(2,l2)))
bBl = VectorOfBlockArrayCoo(blocks,blockids,ax)

dl = lazy_map(*,aBl,bBl)
test_array(dl,[ a*b for (a,b) in zip(aBl,bBl) ])

dl = lazy_map(MulAddMap(2,3),aBl,bBl,bBl)
test_array(dl,[ 2*a*b + 3*c for (a,b,c) in zip(aBl,bBl,bBl) ])

dl = lazy_map(*,aBl,lazy_map(transpose,aBl))
test_array(dl,[ a*transpose(a) for a in aBl ])

end # module
