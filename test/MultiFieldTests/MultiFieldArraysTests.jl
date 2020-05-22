module MultiFieldArraysTests

using Test
using Gridap.Arrays
using Gridap.MultiField
using LinearAlgebra

using Gridap.MultiField: _resize_for_mul!
using Gridap.MultiField: _move_cached_arrays!
using Gridap.MultiField: CachedMultiFieldArray

B1 = 10*ones(Int,3)
B2 = 20*ones(Int,5)
blocks = [B1,B2]
coordinates = [(1,),(2,)]
b = MultiFieldArray(blocks,coordinates)

C1 = 1*ones(Int,3)
C2 = 2*ones(Int,5)
blocks = [C2]
coordinates = [(2,)]
c = MultiFieldArray(blocks,coordinates)

add_to_array!(b,c)
@test b.blocks == [[10, 10, 10], [22, 22, 22, 22, 22]]

A11 = ones(Int,2,3)
A21 = 2*ones(Int,4,3)
A12 = 3*ones(Int,2,5)
blocks = [A11,A21,A12]
coordinates = [(1,1),(2,1),(1,2)]
a = MultiFieldArray(blocks,coordinates)

@test a.blocks == blocks
@test a.coordinates == coordinates
@test a.ptrs == [1 3; 2 0]
@test num_blocks(a) == 4
@test num_stored_blocks(a) == 3
@test get_block_size(a) == (2,2)
@test has_all_blocks(a) == false
@test a[1,1] == A11
@test a[2,1] == A21

B1 = 10*ones(Int,3)
B2 = 20*ones(Int,5)
blocks = [B1,B2]
coordinates = [(1,),(2,)]
b = MultiFieldArray(blocks,coordinates)

c = a*b

@test c.blocks[1] == A11*B1 + A12*B2
@test c.blocks[2] == A21*B1
@test c.coordinates == [(1,),(2,)]

mul!(c,a,b)

@test c.blocks[1] == A11*B1 + A12*B2
@test c.blocks[2] == A21*B1
@test c.coordinates == [(1,),(2,)]

c = a*b
r = CachedMultiFieldArray(c)

_resize_for_mul!(r,a,b)
_move_cached_arrays!(c,r)
mul!(c,10*a,b)

_resize_for_mul!(r,a,b)
_move_cached_arrays!(c,r)
mul!(c,a,b)

A11 = ones(Int,7,3)
A21 = 2*ones(Int,4,3)
A12 = 3*ones(Int,7,5)
blocks = [A11,A21,A12]
coordinates = [(1,1),(2,1),(1,2)]
_a = MultiFieldArray(blocks,coordinates)

_resize_for_mul!(r,_a,b)
_move_cached_arrays!(c,r)
mul!(c,_a,b)

_resize_for_mul!(r,a,b)
_move_cached_arrays!(c,r)
mul!(c,a,b)

add_to_array!(c,c)
add_to_array!(c,10)

end # module
