module BlockFieldArrays


using Gridap.Fields
using Gridap.Arrays
using Gridap.TensorValues
using BlockArrays
using FillArrays
#using Gridap.Fields: MaskedField, MaskedFieldArray

using Test

p = Point(1.0,2.0)
np = 4
x = fill(p,np)
z = fill(p,0)

v1 = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
f1 = MockFieldArray(v1)
v2 = Float64[1,2,3,4]
f2 = MockFieldArray(v2)

axs = (append_ranges([axes(f1)[1],axes(f2)[1]]),)
b1 = BlockFieldArrayCoo(axs,[(1,)],f1)
b2 = BlockFieldArrayCoo(axs,[(2,)],f2)

b1p = evaluate(b1,p)
r = BlockArrayCoo(axs,[(1,)],[v1,])
test_array(b1p,r)

b2p = evaluate(b2,p)
r = BlockArrayCoo(axs,[(2,)],[v2,])
test_array(b2p,r)

b1x = evaluate(b1,x)
test_array(b1x,collect(b1x))

b2x = evaluate(b2,x)
test_array(b2x,collect(b2x))

#using BenchmarkTools
##blockids = [(2,)]
##@btime BlockFieldArrayCoo($axs,$blockids,$f2)
##cache = return_cache(b1,p)
##@btime evaluate!($cache,$b1,$p)
#cache = return_cache(b1,x)
#@btime evaluate!($cache,$b1,$x)

∇b2 = Broadcasting(∇)(b2)
∇b2x = evaluate(∇b2,x)
@test ∇b2x.blockids == b2x.blockids

b2t = transpose(b2)
@test b2t.blockids[1][2] == b2.blockids[1][1]
@test size(b2t) == (1,length(b2))

b2tp = evaluate(b2t,p)
test_array(b2tp,collect(b2tp))

b2tx = evaluate(b2t,x)
test_array(b2tx,collect(b2tx))

∇b2t = Broadcasting(∇)(b2t)
∇b2tx = evaluate(∇b2t,x)
@test ∇b2tx.blockids == b2tx.blockids

# Block of Blocks

axS = (append_ranges([axs[1],axs[1]]),)
b1P = BlockFieldArrayCoo(axS,[(1,)],b1)
b1M = BlockFieldArrayCoo(axS,[(2,)],b1)
b2P = BlockFieldArrayCoo(axS,[(1,)],b2)
b2M = BlockFieldArrayCoo(axS,[(2,)],b2)

b2Pp = evaluate(b2P,p)
@test is_zero_block(b2Pp,2)
@test is_nonzero_block(b2Pp,1)
@test is_zero_block(b2Pp[Block(1)],1)
@test is_nonzero_block(b2Pp[Block(1)],2)
test_array(b2Pp[Block(1)],b2p)

b2Px = evaluate(b2P,x)
@test is_zero_block(b2Px,1,2)
@test is_nonzero_block(b2Px,1,1)
@test is_zero_block(b2Px[Block(1,1)],1,1)
@test is_nonzero_block(b2Px[Block(1,1)],1,2)
test_array(b2Px[Block(1,1)],b2x)

∇b2P = Broadcasting(∇)(b2P)
∇b2Px = evaluate(∇b2P,x)
@test ∇b2Px.blockids == b2Px.blockids

b2Pt = transpose(b2P)
b2Ptp = evaluate(b2Pt,p)

# Global optimizations

ncells = 10
cell_p = Fill(p,ncells)
cell_x = Fill(x,ncells)

cell_f1 = Fill(f1,ncells)
cell_f2 = fill(f2,ncells)

cell_axs = Fill(axs,ncells)

cell_b1 = lazy_map(BlockFieldArrayCooMap(),cell_axs,Fill([(1,)],ncells),cell_f1)
cell_b2 = lazy_map(BlockFieldArrayCooMap(),cell_axs,Fill([(2,)],ncells),cell_f2)

cell_b1p = lazy_map(evaluate,cell_b1,cell_p)
cell_b2p = lazy_map(evaluate,cell_b2,cell_p)

print_op_tree(cell_b2p)

display(cell_b1p[1])
display(cell_b2p[1])


#b1 = BlockFieldArrayCoo(axs,[(1,)],f1)
#b2 = BlockFieldArrayCoo(axs,[(2,)],f2)



#c = Broadcasting(Operation(*))(b2,b2t)
#@show typeof(c)
#display(evaluate(c,p))
#kk



#v = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
#f = MockField.(v)
#
##fi = first(f)
#
#mi = MaskedField(fi,false)
#test_field(mi,p,fi(p))
#test_field(mi,x,fi(x))
#
#mi = MaskedField(fi,true)
#test_field(mi,p,0*fi(p))
#test_field(mi,x,0*fi(x))
#
#m = MaskedFieldArray(axes(f),f,false)
#test_field_array(m,p,evaluate(f,p))
#test_field_array(m,x,evaluate(f,x))
#
#m = MaskedFieldArray(axes(f),f,true)
#test_field_array(m,p,0*evaluate(f,p))
#test_field_array(m,x,0*evaluate(f,x))
#
#
#v1 = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
#f1 = MockField.(v1)
#v2 = Float64[1,2,3,4]
#f2 = MockField.(v2)
#
#axs = (append_ranges([axes(f1)[1],axes(f2)[1]]),)
#
#f11 = MaskedFieldArray(axes(f1),f1,false)
#f12 = MaskedFieldArray(axes(f2),f1,true)
#
#b1 = BlockArrayCoo(axs,[(1,)],[f11],[1,-1],[f12])
#
#display(b1)
#
#display(evaluate(b1,p))
#
#
#kk



#display(m)
##test_field_array(m,p,evaluate(f,p))
##test_field_array(m,x,evaluate(f,x))
#
#m = MaskedFieldArray(axes(f),f,true)
##test_field_array(m,p,0*evaluate(f,p))
##test_field_array(m,x,0*evaluate(f,x))




#v1 = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
#f1 = MockField.(v1)
#
#v2 = Float64[1,2,3]
#f2 = MockField.(v1)


end # module
