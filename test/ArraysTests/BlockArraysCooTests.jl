module BlockArraysCooTests

using Test
using Gridap.Arrays
using Gridap.Algebra
using BlockArrays
using LinearAlgebra
using Gridap.Arrays: BlockArrayCooMap

# Check MultiLevelBlockedUnitRange

r1 = Base.OneTo(2)
r2 = Base.OneTo(3)

r = append_ranges([r1,r2])
@test length(r) == 5
@test num_blocks_equal(r,r)

s1 = Base.OneTo(1)
s2 = Base.OneTo(2)
s3 = Base.OneTo(0)

s = append_ranges([s1,s2,s3])
@test length(s) == 3
@test !num_blocks_equal(r,s)

o = append_ranges([r,s])
@test length(o) == 8
@test blocks_equal(o,o)
@test num_blocks_equal(o,o)
@test o.local_ranges[1] === r
@test o.local_ranges[2] === s

a = append_ranges([r,s])
@test blocks_equal(o,a)
@test num_blocks_equal(o,a)

#using BenchmarkTools
#@btime blocks_equal($o,$a)
#@btime num_blocks_equal($o,$a)

b = append_ranges([s,s])
@test !blocks_equal(a,b)
@test BlockArrays.blocksize(a) == BlockArrays.blocksize(b) # This only refers to the top level blocks
@test !num_blocks_equal(a,b) # This takes into account all levels

r1_axes = (r1,)
r2_axes = (r2,)
s1_axes = (s1,)
s2_axes = (s2,)

r_axes = (append_ranges(map(first,[r1_axes,r2_axes])),)
s_axes = (append_ranges(map(first,[s1_axes,s2_axes])),)
o_axes = (append_ranges(map(first,[r_axes,s_axes])),)
@test blocks_equal(o_axes,o_axes)
@test num_blocks_equal(o_axes,o_axes)

#using BenchmarkTools
#@btime blocks_equal($o_axes,$o_axes)
#@btime num_blocks_equal($o_axes,$o_axes)

# BlockArrayCoo

# Minimal constructor

blocks = Matrix{Float64}[ [1 2; 3 4], [5 6; 7 8; 9 10] ]
blockids = [(1,1),(2,1)]
ax = (blockedrange([2,3]), blockedrange([2,4]))
a = BlockArrayCoo(ax,blockids,blocks)
a = BlockArrayCooMap((2,2),blockids)(ax,blocks...)

k = BlockArrayCooMap((2,2),blockids)
@test blocksize(k) == (2,2)
@test is_zero_block(k,1,1) == false
@test is_zero_block(k,2,1) == false
@test is_zero_block(k,1,2) == true
@test is_zero_block(k,2,2) == true
@test is_nonzero_block(k,1,1) == true
@test is_nonzero_block(k,first(eachblockid(k)))

# Specific API

@test is_zero_block(rand(4,5),1,1) == false
@test is_zero_block(a,2,2) == true
@test is_zero_block(a,2,1) == false
@test is_nonzero_block(a,2,1) == true
@test is_zero_block(a,Block(1,2)) == true
@test is_zero_block(a,Block(1),Block(2)) == true

# AbstractBlockArray interface

@test a[Block(1),Block(1)] === blocks[1]
@test a[Block(1,1)] === blocks[1]
@test a[Block(1,2)] === a.zero_blocks[1]

c = copy(a[Block(1,1)])
fill!(c,0)

getblock!(c,a,1,1)
@test c == a[Block(1,1)]
@test a[Block(1,1)] === blocks[1]

c = copy(a)
@test isa(c,BlockArrayCoo)
@test c == a
c[Block(1,1)] = a[Block(1,1)]
@test c[Block(1,1)] === a[Block(1,1)]
@test a[Block(1,1)] === blocks[1]
c[1,2] = 3
@test c[1,2] == 3

#using BenchmarkTools
#@btime getindex($c,1,2)
#bi = BlockIndex((1,1),(1,2))
#@btime getindex($c,$bi)

@test a[BlockIndex((1,1),(2,1))] === blocks[1][2,1]
@test a[BlockIndex(1,2),BlockIndex(1,1)] === blocks[1][2,1]
@test a[2,1] === blocks[1][2,1]
@test a[3,2] === blocks[2][1,2]
@test axes(a) === ax
@test size(a) == (5,6)

for (i,b) in enumerateblocks(a)
  @test a[i] === b
end

b21 = zeros(3,2)
getblock!(b21,a,Block(2,1))
@test b21 == a[Block(2,1)]

b12 = ones(2,4)
getblock!(b12,a,Block(1,2))
@test b12 == a[Block(1,2)]

blocks = Vector{Float64}[ [1,2,3] ]
blockids = [(2,)]
ax = (blockedrange([2,3,4]),)
a = BlockArrayCoo(ax,blockids,blocks)

@test a[Block(1)] === a.zero_blocks[1]
@test a[Block(2)] === blocks[1]
@test a[Block(3)] === a.zero_blocks[2]
@test a[BlockIndex(2,3)] === blocks[1][3]

blocks = Matrix{Float64}[ [1 2; 3 4], [5 6 7 8; 9 10 11 12; 13 14 15 16], [1 2 3 4; 5 6 7 8], [1 2 3; 4 5 6; 7 8 9] ]
blockids = [(1,1),(2,2),(1,2),(3,3)]
ax = (blockedrange([2,3,3]), blockedrange([2,4,3]))
a = BlockArrayCoo(ax,blockids,blocks)

blocks = Vector{Float64}[ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
ax = (blockedrange([2,4,3]),)
b = BlockArrayCoo(ax,blockids,blocks)

# similar

c = similar(a)
@test typeof(c) == typeof(a)
@test axes(a) == axes(c)
@test a.blockids == c.blockids

c = similar(a,Int)
@test axes(a) == axes(c)
@test a.blockids == c.blockids

ax = (blockedrange([2,3,3]), blockedrange([2,4,3]))
c = similar(a,eltype(a),ax)
@test typeof(c) == typeof(a)
@test axes(a) == axes(c)
@test a.blockids == c.blockids

# In this case, we can still preserve the zero block structure
ax = (blockedrange([2,5,3]), blockedrange([2,4,1]))
c = similar(a,eltype(a),ax)
@test typeof(c) == typeof(a)
@test axes(a) != axes(c)
@test axes(c) === ax
@test a.blockids == c.blockids
fill_entries!(c,2)

# In this case, we cannot preserve the zero block structure
ax = (blockedrange([2,5,3,5]), blockedrange([2,4,1]))
c = similar(a,eltype(a),ax)
@test typeof(c) == typeof(a)
@test axes(a) != axes(c)
@test axes(c) === ax
@test a.blockids != c.blockids
@test length(c.blockids) == 12

# zero

c = similar(a)
fill_entries!(c,1)

#using BenchmarkTools
#@btime fill_entries!($c,1)

z = zero(c)
@test isa(z,BlockArrayCoo)
@test length(z.blocks) == 0

# + and -

c = z + a
@test isa(c,BlockArrayCoo)
@test c == a
@test Array(c) == Array(a)

c = z - a
@test isa(c,BlockArrayCoo)
@test c == -a
@test Array(c) == Array(-a)

# mat vec *

c = a*b
@test isa(c,BlockArrayCoo)
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

#using BenchmarkTools
#@btime mul!($c,$a,$b)
#@btime mul!($d,$a,$b,2,3)

# mat mat *

b = transpose(a)
@test isa(b,BlockArrayCoo)
c = a*b
@test isa(c,BlockArrayCoo)
@test axes(c,1) === axes(a,1)
@test axes(c,2) === axes(a,1)
@test blocksize(c) == (3,3)
@test Array(a)*Array(b) == c

mul!(c,a,b)
@test axes(c,1) === axes(a,1)
@test axes(c,2) === axes(a,1)
@test blocksize(c) == (3,3)
@test Array(a)*Array(b) == c

#using BenchmarkTools
#@btime mul!($c,$a,$b)

# Cached array
cc = CachedArray(c)

# Compatible axes, we do not modify the cache
axs = (blockedrange([2,3,3]), blockedrange([2,3,3]))
@test blocks_equal(axs,axes(c))
setaxes!(cc,axs)
@test map(blocklasts,axes(cc.array)) == map(blocklasts,axs)
@test cc.array === c

#using BenchmarkTools
#@btime setaxes!($cc,$axs)

# Not compatible axes but same number of blocks
# we create a new cache, but we can preserve the zero block structure
axs = (blockedrange([4,5,3]), blockedrange([2,4,3]))
@test !blocks_equal(axs,axes(c))
@test num_blocks_equal(axs,axes(c))
setaxes!(cc,axs)
@test map(blocklasts,axes(cc.array)) == map(blocklasts,axs)
@test cc.array != c
@test cc.array.blockids == c.blockids
scale_entries!(cc.array,2)

# Not compatible axes and different number of blocks
# we create a new cache without preserving zero block structure
axs = (blockedrange([5,4,3,4]), blockedrange([4,2,3]))
@test !blocks_equal(axs,axes(cc.array))
@test !num_blocks_equal(axs,axes(cc.array))
setaxes!(cc,axs)
@test map(blocklasts,axes(cc.array)) == map(blocklasts,axs)
@test cc.array != c
@test cc.array.blockids != c.blockids

axs1 = (blockedrange([5,4,3]), blockedrange([4,2,3]))
axs2 = (blockedrange([4,5,3]), blockedrange([2,4,3]))
axs3 = (blockedrange([4,5,3]), blockedrange([2,4,3]))
@test blocks_equal(axs1,axs2) == false
@test blocks_equal(axs2,axs2)
@test blocks_equal(axs2,axs3)

blocks = [ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
axs = (blockedrange([2,4,3]),)
b = BlockArrayCoo(axs,blockids,blocks)

cb = CachedArray(b)

setaxes!(cb,axs)
@test map(blocklasts,axes(cb.array)) == map(blocklasts,axs)
@test cb.array === b

axs = (blockedrange([3,2,3]),)
setaxes!(cb,axs)
@test map(blocklasts,axes(cb.array)) == map(blocklasts,axs)
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
a = BlockArrayCoo(ax,blockids,blocks)

blocks = [ [1 2; 3 4], [1 2; 4 5; 8 9], [5 6 7 8; 9 10 11 12; 13 14 15 16] ]
blockids = [(1,1),(2,1),(3,2)]
a = BlockArrayCoo(ax,blockids,blocks)

c = a - 2*a
@test isa(c,BlockArrayCoo)
@test Array(a)-2*Array(a) == c

@test a[1,2,1] == a[1,2]

blocks = [ 10*[1,2], 20*[1,2,3] ]
blockids = [(1,),(3,)]
axs = (blockedrange([2,4,3]),)
b = BlockArrayCoo(axs,blockids,blocks)

@test b[2] == b[2,1]

# Blocks of Blocks

axs = (blockedrange([2,4,3]),)
bL = BlockArrayCoo(axs,[(1,)],[[1,2]])
bR = BlockArrayCoo(axs,[(1,),(3,)],[[2,3],[4,5,6]])

axsS = (append_ranges(map(first,[axs,axs])),)

bS = BlockArrayCoo(axsS,[(1,),(2,)],[bL,bR])

zS = zero(bS)
@test isa(zS,BlockArrayCoo)
@test isa(zS[Block(1)],BlockArrayCoo)
@test isa(zS[Block(2)],BlockArrayCoo)
@test length(zS.blocks) == 0
@test length(zS[Block(1)].blocks) == 0
@test length(zS[Block(2)].blocks) == 0

cS = 2*bS
@test isa(cS,BlockArrayCoo)
@test isa(cS[Block(1)],BlockArrayCoo)

dS = BlockArrayCoo(axsS,[(2,)],[bR])

r = cS - dS
@test isa(r,BlockArrayCoo)
@test isa(r[Block(1)],BlockArrayCoo)

r = cS + dS
@test isa(r,BlockArrayCoo)
@test isa(r[Block(1)],BlockArrayCoo)

# Block of block matrices

ax = (blockedrange([2,3]), blockedrange([2,4]))
aLL = BlockArrayCoo(ax,[(1,1),(2,1)],[[1 2; 3 4],[5 6; 7 8; 9 10]])
aLR = BlockArrayCoo(ax,[(2,2),(2,1)],[[1 2 5 6; 3 4 1 2; 1 2 3 4],[5 6; 7 8; 9 10]])
aRR = BlockArrayCoo(ax,[(1,1),(2,1)],[[1 2; 3 4],[5 6; 7 8; 9 10]])

axsS = ( append_ranges([ax[1],ax[1]]) , append_ranges([ax[2],ax[2]]))
aS = BlockArrayCoo(axsS,[(1,1),(1,2),(2,2)],[aLL,aLR,aRR])
@test is_zero_block(aS,2,1)
@test is_nonzero_block(aS,2,2)

aS2 = copy(aS)
@test isa(aS2,BlockArrayCoo)
@test isa(aS2[Block(2),Block(2)],BlockArrayCoo)
@test is_zero_block(aS2,2,1)
@test is_nonzero_block(aS2,2,2)

aSt = transpose(aS)
@test Array(aSt) == transpose(Array(aS))
@test isa(aSt,BlockArrayCoo)
@test isa(aSt[Block(2),Block(2)],BlockArrayCoo)
@test is_zero_block(aSt,1,2)
@test is_nonzero_block(aSt,2,2)

rS = aS*aSt
@test rS == Array(aS)*Array(aSt)
@test isa(rS,BlockArrayCoo)
@test isa(rS[Block(2),Block(2)],BlockArrayCoo)

mul!(rS,aS,aSt)
@test rS == Array(aS)*Array(aSt)
@test isa(rS,BlockArrayCoo)
@test isa(rS[Block(2),Block(2)],BlockArrayCoo)

#using BenchmarkTools
#@btime mul!($rS,$aS,$aSt)

axs =(blockedrange([2,4]),)
bL = BlockArrayCoo(axs,[(1,)],[[1,2]])
bR = BlockArrayCoo(axs,[(1,),(2,)],[[2,3],[4,5,6,8]])
axsS = (append_ranges(map(first,[axs,axs])),)
bS = BlockArrayCoo(axsS,[(1,),(2,)],[bL,bR])

rS = aS*bS
@test rS == Array(aS)*Array(bS)
@test isa(rS,BlockArrayCoo)
@test isa(rS[Block(2)],BlockArrayCoo)

mul!(rS,aS,bS)
@test rS == Array(aS)*Array(bS)
@test isa(rS,BlockArrayCoo)
@test isa(rS[Block(2)],BlockArrayCoo)

#using BenchmarkTools
#@btime mul!($rS,$aS,$bS)

cS = copy(rS)
mul!(rS,aS,bS,3,2)
@test rS == 3*Array(aS)*Array(bS) + 2*cS
@test isa(rS,BlockArrayCoo)
@test isa(rS[Block(2)],BlockArrayCoo)

#using BenchmarkTools
#@btime mul!($rS,$aS,$bS,3,2)

#using BenchmarkTools
#@btime Arrays._same_axes($axsA,$axsA)

end # module
