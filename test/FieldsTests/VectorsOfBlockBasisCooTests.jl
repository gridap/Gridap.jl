module InsertInBlocksTests

using Test
using BlockArrays
using FillArrays
using Gridap.Arrays
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Fields: MockField, MockBasis

l = 10
np = 3
ndofs1 = 4
ndofs2 = 6

p = Point(1,2)
x = fill(p,np)
z = 2.0

v = VectorValue(3.0,1.5)
w = VectorValue(3.4,3.5)
a = MockBasis{2}(v,ndofs2)
b = MockBasis{2}(w,ndofs1)

al = Fill(a,l)
atl = trialize_array_of_bases(al)
bl = Fill(b,l)
xl = Fill(x,l)

ax1 = Fill((Base.OneTo(ndofs1),),l)
ax2 = Fill((Base.OneTo(ndofs2),),l)
aBl = insert_array_of_bases_in_block(al,ax1,ax2,2)
aBl_x = evaluate(aBl,xl)
@test isa(aBl,VectorOfBlockBasisCoo)
@test isa(aBl_x,VectorOfBlockArrayCoo)

bBl = insert_array_of_bases_in_block(bl,ax1,ax2,1)
bBl_x = evaluate(bBl,xl)
@test isa(bBl,VectorOfBlockBasisCoo)
@test isa(bBl_x,VectorOfBlockArrayCoo)

ax1 = Fill((Base.OneTo(1),Base.OneTo(ndofs1)),l)
ax2 = Fill((Base.OneTo(1),Base.OneTo(ndofs2)),l)
atBl = insert_array_of_bases_in_block(atl,ax1,ax2,2)
atBl_x = evaluate(atBl,xl)
@test isa(atBl,VectorOfBlockBasisCoo)
@test isa(atBl_x,VectorOfBlockArrayCoo)

ax1 = aBl.axes
ax2 = aBl.axes
aSl = insert_array_of_bases_in_block(aBl,ax1,ax2,1)
@test isa(aSl,VectorOfBlockBasisCoo)
aSl_x = evaluate(aSl,xl)
@test isa(aSl_x,VectorOfBlockArrayCoo)

ax1 = atBl.axes
ax2 = atBl.axes
atSl = insert_array_of_bases_in_block(atBl,ax1,ax2,1)
@test isa(atSl,VectorOfBlockBasisCoo)
atSl_x = evaluate(atSl,xl)
@test isa(atSl_x,VectorOfBlockArrayCoo)

blocks = (al,)
blockids = [(2,)]
ranges_i = (blockedrange([ndofs1,ndofs2]),)
ranges = Fill(ranges_i,l)
aBl = VectorOfBlockBasisCoo(blocks,blockids,ranges)

blocks = (bl,)
blockids = [(1,)]
bBl = VectorOfBlockBasisCoo(blocks,blockids,ranges)

blocks = (atl,)
blockids = [(1,2)]
ranges_i = (blockedrange([1]),blockedrange([ndofs1,ndofs2]))
ranges = Fill(ranges_i,l)
atBl = VectorOfBlockBasisCoo(blocks,blockids,ranges)

aBl_x = evaluate(aBl,xl)
@test isa(aBl_x,VectorOfBlockArrayCoo)
@test isa(aBl_x.axes,Fill)
@test aBl_x.axes[1] == (blockedrange([np]),blockedrange([ndofs1,ndofs2]))

bBl_x = evaluate(bBl,xl)
@test isa(bBl_x,VectorOfBlockArrayCoo)

atBl_x = evaluate(atBl,xl)
@test isa(atBl_x,VectorOfBlockArrayCoo)
@test isa(atBl_x.axes,Fill)
@test atBl_x.axes[1] == (blockedrange([np]),blockedrange([1]),blockedrange([ndofs1,ndofs2]))

atBl = trialize_array_of_bases(aBl)
atBl_x = evaluate(atBl,xl)
@test isa(atBl_x,VectorOfBlockArrayCoo)

btBl = trialize_array_of_bases(bBl)
btBl_x = evaluate(btBl,xl)
@test isa(btBl_x,VectorOfBlockArrayCoo)

∇aBl = ∇(aBl)
@test isa(∇aBl,VectorOfBlockBasisCoo)
∇aBl_x = evaluate(∇aBl,xl)
@test isa(∇aBl_x,VectorOfBlockArrayCoo)

cl = operate_arrays_of_fields(-,aBl,bBl)
cl_x = evaluate(cl,xl)
@test isa(cl_x,VectorOfBlockArrayCoo)
∇cl = ∇(cl)
∇cl_x = evaluate(∇cl,xl)
@test isa(cl_x,VectorOfBlockArrayCoo)

cl = operate_arrays_of_fields(⋅,atBl,cl)
cl_x = evaluate(cl,xl)
@test isa(cl_x,VectorOfBlockArrayCoo)
@test is_zero_block(cl_x,Block(1,1,1))
@test is_nonzero_block(cl_x,Block(1,1,2))
@test is_zero_block(cl_x,Block(1,2,1))
@test is_nonzero_block(cl_x,Block(1,2,2))

ids = [2,3,5,1,2]
cj = reindex(aBl,ids)
xj = reindex(xl,ids)
@test isa(cj,VectorOfBlockBasisCoo)
cj_x = evaluate(cj,xj)
@test isa(cj_x.axes,Fill)
@test cj_x.axes[1] == (blockedrange([np]),blockedrange([ndofs1,ndofs2]))
@test length(cj_x) == length(ids)
@test length(cj) == length(ids)

v = VectorValue(1.0,2.0,3.0)
f = MockField{2}(v)
fl = Fill(f,l)

cl = compose_field_arrays(aBl,fl)
@test isa(cl,VectorOfBlockBasisCoo)


#cache = array_cache(cl_x)
#using BenchmarkTools
#@btime getindex!($cache,$cl_x,3)

end # module
