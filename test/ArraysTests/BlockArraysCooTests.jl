module BlockArraysCooTests

using Test
using Gridap.Arrays
using BlockArrays
using LinearAlgebra

blocks = [ [1 2; 3 4], [5 6; 7 8; 9 10] ]
blockids = [(1,1),(2,1)]
ax = (blockedrange([2,3]), blockedrange([2,4]))

a = BlockArrayCoo(blocks,blockids,ax)
@test a[Block(1),Block(1)] === blocks[1]
@test a[Block(1,1)] === blocks[1]
@test a[Block(1,2)] === a.zero_blocks[1]
@test a[BlockIndex((1,1),(2,1))] === blocks[1][2,1]
@test a[BlockIndex(1,2),BlockIndex(1,1)] === blocks[1][2,1]
@test a[2,1] === blocks[1][2,1]
@test a[3,2] === blocks[2][1,2]
@test axes(a) === ax
@test size(a) == (5,6)

@test is_zero_block(a,2,2) == true
@test is_zero_block(a,2,1) == false
@test is_nonzero_block(a,2,1) == true
@test is_zero_block(a,Block(1,2)) == true
@test is_zero_block(a,Block(1),Block(2)) == true

for (i,b) in enumerateblocks(a)
  @test a[i] === b
end

b21 = zeros(3,2)
getblock!(b21,a,Block(2,1))
@test b21 == a[Block(2,1)]

b12 = ones(2,4)
getblock!(b12,a,Block(1,2))
@test b12 == a[Block(1,2)]

blocks = [ [1,2,3] ]
blockids = [(2,)]
ax = (blockedrange([2,3,4]),)
a = BlockArrayCoo(blocks,blockids,ax)

@test a[Block(1)] === a.zero_blocks[1]
@test a[Block(2)] === blocks[1]
@test a[Block(3)] === a.zero_blocks[2]
@test a[BlockIndex(2,3)] === blocks[1][3]

blocks = [ [1 2; 3 4], [5 6 7 8; 9 10 11 12; 13 14 15 16], [1 2 3 4; 5 6 7 8], [1 2 3; 4 5 6; 7 8 9] ]
blockids = [(1,1),(2,2),(1,2),(3,3)]
ax = (blockedrange([2,3,3]), blockedrange([2,4,3]))
a = BlockArrayCoo(blocks,blockids,ax)

blocks = [ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
ax = (blockedrange([2,4,3]),)
b = BlockArrayCoo(blocks,blockids,ax)

c = a*b
@test axes(c,1) === axes(a,1)
@test blocksize(c) == (3,)
@test Array(a)*Array(b) == c

mul!(c,a,b)
@test axes(c,1) === axes(a,1)
@test blocksize(c) == (3,)
@test Array(a)*Array(b) == c

d = copy(c)
mul!(d,a,b,2,3)
@test axes(d,1) === axes(a,1)
@test blocksize(d) == (3,)
@test 2*Array(a)*Array(b) + 3*Array(c) == d

b = Transpose(a)
c = a*b
@test axes(c,1) === axes(a,1)
@test axes(c,2) === axes(a,1)
@test blocksize(c) == (3,3)
@test Array(a)*Array(b) == c

mul!(c,a,b)
@test axes(c,1) === axes(a,1)
@test axes(c,2) === axes(a,1)
@test blocksize(c) == (3,3)
@test Array(a)*Array(b) == c

cc = CachedArray(c)

axs = (blockedrange([2,3,3]), blockedrange([2,3,3]))
setaxes!(cc,axs)
@test cc.array === c

axs = (blockedrange([4,5,3]), blockedrange([2,4,3]))
setaxes!(cc,axs)
fill!(cc.array,0)

blocks = [ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
axs = (blockedrange([2,4,3]),)
b = BlockArrayCoo(blocks,blockids,axs)

cb = CachedArray(b)

setaxes!(cb,axs)
@test cb.array === b

axs = (blockedrange([3,2,3]),)
setaxes!(cb,axs)
@test size(cb) == (8,)

c = copy(a)
@test isa(c,BlockArrayCoo)
@test c == a
@test axes(c) == axes(a)

fill!(c,0)
copyto!(c,a)
@test c == a
@test axes(c) == axes(a)

c = 2*a
@test isa(c,BlockArrayCoo)
@test 2*Array(a) == c

c = a*2
@test isa(c,BlockArrayCoo)
@test 2*Array(a) == c

c = a + a
@test isa(c,BlockArrayCoo)
@test 2*Array(a) == c

blocks = [ [1 2; 3 4], [5 6 7 8; 9 10 11 12; 13 14 15 16], [1 2 3 4; 5 6 7 8], [1 2 3; 4 5 6; 7 8 9] ]
blockids = [(1,1),(2,2),(1,2),(3,3)]
ax = (blockedrange([2,3,3]), blockedrange([2,4,3]))
a = BlockArrayCoo(blocks,blockids,ax)

blocks = [ [1 2; 3 4], [1 2; 4 5; 8 9], [5 6 7 8; 9 10 11 12; 13 14 15 16] ]
blockids = [(1,1),(2,1),(3,2)]
b = BlockArrayCoo(blocks,blockids,ax)

c = a + b
@test isa(c,BlockArrayCoo)
@test Array(a)+Array(b) == c

c = a - 2*b
@test isa(c,BlockArrayCoo)
@test Array(a)-2*Array(b) == c

@test a[1,2,1] == a[1,2]

blocks = [ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
axs = (blockedrange([2,4,3]),)
b = BlockArrayCoo(blocks,blockids,axs)

@test b[2] == b[2,1]

end # module
