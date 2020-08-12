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
bl = Fill(b,l)
xl = Fill(x,l)

blocks = (al,)
blockids = [(1,2)]
axs_i = (blockedrange([np]),blockedrange([ndofs1,ndofs2]))
axs = Fill(axs_i,l)
aBl = VectorOfBlockBasisCoo(blocks,blockids,axs)

blocks = (bl,)
blockids = [(1,1)]
bBl = VectorOfBlockBasisCoo(blocks,blockids,axs)

aBl_x = evaluate(aBl,xl)
@test isa(aBl_x,VectorOfBlockArrayCoo)

bBl_x = evaluate(bBl,xl)
@test isa(bBl_x,VectorOfBlockArrayCoo)

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

#cache = array_cache(cl_x)
#using BenchmarkTools
#@btime getindex!($cache,$cl_x,3)

end # module
