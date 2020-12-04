module BlockFieldArrays


using Gridap.Fields
using Gridap.Arrays
using Gridap.TensorValues
using BlockArrays
using FillArrays
#using Gridap.Fields: MaskedField, MaskedFieldArray

using Gridap.Fields: BlockFieldArrayCoo
using Gridap.Fields: BlockFieldArrayCooMap
using Gridap.Arrays: BlockArrayCooMap

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

cell_b1 = lazy_map(BlockFieldArrayCooMap((2,),[(1,)]),cell_axs,cell_f1)
cell_b2 = lazy_map(BlockFieldArrayCooMap((2,),[(2,)]),cell_axs,cell_f2)

cell_b1p = lazy_map(evaluate,cell_b1,cell_p)
cell_b2p = lazy_map(evaluate,cell_b2,cell_p)
@test isa(cell_b2p.g.value,BlockArrayCooMap)

cell_b1x = lazy_map(evaluate,cell_b1,cell_x)
cell_b2x = lazy_map(evaluate,cell_b2,cell_x)
@test isa(cell_b2x.g.value,BlockArrayCooMap)

cell_∇b2 = lazy_map(Broadcasting(∇),cell_b2)
@test isa(cell_∇b2.g.value,BlockFieldArrayCooMap)

phi(x) = 2*x
cell_phi = Fill(GenericField(phi),ncells)
cell_b2ophi = lazy_map(Broadcasting(∘),cell_b2,cell_phi)
@test isa(cell_b2ophi.g.value,BlockFieldArrayCooMap)

face_cell = [3,2,1,2]
face_b2 = lazy_map(Reindex(cell_b2),face_cell)
@test isa(face_b2.g.value,BlockFieldArrayCooMap)

cell_b2t = lazy_map(transpose,cell_b2)
@test isa(cell_b2t.g.value,BlockFieldArrayCooMap)

cell_b2tx = lazy_map(evaluate,cell_b2t,cell_x)
test_array(cell_b2tx,Fill(b2tx,ncells))

cell_g2 = lazy_map(Broadcasting(Operation(-)),cell_b2)
cell_g2x = lazy_map(evaluate,cell_g2,cell_x)
@test isa(cell_g2x.g.value,BlockArrayCooMap)
test_array(cell_g2x,-collect(cell_b2x))

cell_f = Fill(ConstantField(2.0),ncells)
cell_g2 = lazy_map(Broadcasting(Operation(*)),cell_b2,cell_f)
cell_g2x = lazy_map(evaluate,cell_g2,cell_x)
@test isa(cell_g2x.g.value,BlockArrayCooMap)
r = map(i->2.0*i,cell_b2x)
test_array(cell_g2x,r)

cell_f = Fill(ConstantField(VectorValue(1,1)),ncells)
cell_g1 = lazy_map(Broadcasting(Operation(⋅)),cell_b1,cell_f)
cell_g1x = lazy_map(evaluate,cell_g1,cell_x)
@test isa(cell_g1x.g.value,BlockArrayCooMap)

cell_g = lazy_map(Broadcasting(Operation(-)),cell_g1,cell_g2)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)

cell_g = lazy_map(Broadcasting(Operation(*)),cell_g1,cell_b2t)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)

cell_g = lazy_map(Broadcasting(Operation(*)),cell_b2t,cell_g1)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)

ϕ = GenericField(x->2*x)
cell_ϕ = fill(ϕ,ncells)
cell_J = lazy_map(∇,cell_ϕ)
w = rand(np)
cell_w = Fill(w,ncells)

cell_i = lazy_map(integrate,cell_b2,cell_x,cell_w,cell_J)
@test isa(cell_i.g.value,BlockArrayCooMap)

cell_g = lazy_map(Broadcasting(Operation(*)),cell_g1,cell_b2t)
cell_i = lazy_map(integrate,cell_g,cell_x,cell_w,cell_J)
@test isa(cell_i.g.value,BlockArrayCooMap)

f(a) = 2*a
cell_g = lazy_map(Broadcasting(Operation(f)),cell_b1)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
test_array(cell_gx,map((a)->f.(a),cell_b1x))

cell_f = Fill(ConstantField(3),ncells)
cell_fx = lazy_map(evaluate,cell_f,cell_x)
f(a,b) = 2*a + a*b
cell_g = lazy_map(Broadcasting(Operation(f)),cell_b2,cell_f)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
test_array(cell_gx,map((a,b)->f.(a,b),cell_b2x,cell_fx))

f(b,a) = 2*a + a*b
cell_g = lazy_map(Broadcasting(Operation(f)),cell_f,cell_b2)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
r = map((a,b)->f.(a,b),cell_fx,cell_b2x)
test_array(cell_gx,r)

# Blocks of blokcs

cell_axsS = lazy_map(cell_axs) do axs
  (append_ranges([axs[1],axs[1]]),)
end

cell_b1L = lazy_map(BlockFieldArrayCooMap((2,),[(1,)]),cell_axsS,cell_b1)
cell_b1R = lazy_map(BlockFieldArrayCooMap((2,),[(2,)]),cell_axsS,cell_b1)
cell_b2L = lazy_map(BlockFieldArrayCooMap((2,),[(1,)]),cell_axsS,cell_b2)
cell_b2R = lazy_map(BlockFieldArrayCooMap((2,),[(2,)]),cell_axsS,cell_b2)

cell_b1Lx = lazy_map(evaluate,cell_b1L,cell_x)
cell_b1Rx = lazy_map(evaluate,cell_b1R,cell_x)
cell_b2Lx = lazy_map(evaluate,cell_b2L,cell_x)
cell_b2Rx = lazy_map(evaluate,cell_b2R,cell_x)
@test isa(cell_b1Lx.g.value,BlockArrayCooMap)
@test isa(cell_b2Rx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_b1Lx[1],1,1) == false
@test is_zero_block(cell_b1Lx[1],1,2) == true
@test is_zero_block(cell_b1Lx[1][Block(1,1)],1,1) == false
@test is_zero_block(cell_b1Lx[1][Block(1,1)],1,2) == true
@test is_zero_block(cell_b2Lx[1],1,1) == false
@test is_zero_block(cell_b2Lx[1],1,2) == true
@test is_zero_block(cell_b2Lx[1][Block(1,1)],1,1) == true
@test is_zero_block(cell_b2Lx[1][Block(1,1)],1,2) == false

cell_∇b1R = lazy_map(Broadcasting(∇),cell_b1R)
cell_∇b2L = lazy_map(Broadcasting(∇),cell_b2L)
cell_∇b1Rx = lazy_map(evaluate,cell_∇b1R,cell_x)
cell_∇b2Lx = lazy_map(evaluate,cell_∇b2L,cell_x)
@test isa(cell_∇b1Rx.g.value,BlockArrayCooMap)
@test isa(cell_∇b2Lx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_∇b1Rx[1],1,1) == true
@test is_zero_block(cell_∇b1Rx[1],1,2) == false
@test is_zero_block(cell_∇b1Rx[1][Block(1,2)],1,1) == false
@test is_zero_block(cell_∇b1Rx[1][Block(1,2)],1,2) == true
@test is_zero_block(cell_∇b2Lx[1],1,1) == false
@test is_zero_block(cell_∇b2Lx[1],1,2) == true
@test is_zero_block(cell_∇b2Lx[1][Block(1,1)],1,1) == true
@test is_zero_block(cell_∇b2Lx[1][Block(1,1)],1,2) == false

cell_b1Rt = lazy_map(transpose,cell_b1R)
cell_b2Lt = lazy_map(transpose,cell_b2L)
cell_b2Rt = lazy_map(transpose,cell_b2R)
cell_b1Rtx = lazy_map(evaluate,cell_b1Rt,cell_x)
cell_b2Ltx = lazy_map(evaluate,cell_b2Lt,cell_x)
cell_b2Rtx = lazy_map(evaluate,cell_b2Rt,cell_x)
@test isa(cell_b1Rtx.g.value,BlockArrayCooMap)
@test isa(cell_b2Ltx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_b1Rtx[1],1,1,1) == true
@test is_zero_block(cell_b1Rtx[1],1,1,2) == false
@test is_zero_block(cell_b1Rtx[1][Block(1,1,2)],1,1,1) == false
@test is_zero_block(cell_b1Rtx[1][Block(1,1,2)],1,1,2) == true
@test is_zero_block(cell_b2Ltx[1],1,1,1) == false
@test is_zero_block(cell_b2Ltx[1],1,1,2) == true
@test is_zero_block(cell_b2Ltx[1][Block(1,1,1)],1,1,1) == true
@test is_zero_block(cell_b2Ltx[1][Block(1,1,1)],1,1,2) == false

cell_∇b1Rt = lazy_map(Broadcasting(∇),cell_b1Rt)
cell_∇b2Lt = lazy_map(Broadcasting(∇),cell_b2Lt)
cell_∇b1Rtx = lazy_map(evaluate,cell_∇b1Rt,cell_x)
cell_∇b2Ltx = lazy_map(evaluate,cell_∇b2Lt,cell_x)
@test isa(cell_∇b1Rtx.g.value,BlockArrayCooMap)
@test isa(cell_∇b2Ltx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_∇b1Rtx[1],1,1,1) == true
@test is_zero_block(cell_∇b1Rtx[1],1,1,2) == false
@test is_zero_block(cell_∇b1Rtx[1][Block(1,1,2)],1,1,1) == false
@test is_zero_block(cell_∇b1Rtx[1][Block(1,1,2)],1,1,2) == true
@test is_zero_block(cell_∇b2Ltx[1],1,1,1) == false
@test is_zero_block(cell_∇b2Ltx[1],1,1,2) == true
@test is_zero_block(cell_∇b2Ltx[1][Block(1,1,1)],1,1,1) == true
@test is_zero_block(cell_∇b2Ltx[1][Block(1,1,1)],1,1,2) == false

f(a) = 3*a
cell_g = lazy_map(Broadcasting(Operation(f)),cell_b2Lt)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_gx[1],1,1,1) == false
@test is_zero_block(cell_gx[1],1,1,2) == true
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,1) == true
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,2) == false
r = map(i->3.0*i,cell_b2Ltx)
test_array(cell_gx,r)

cell_f = Fill(ConstantField(3),ncells)
cell_fx = lazy_map(evaluate,cell_f,cell_x)
cell_g = lazy_map(Broadcasting(Operation(*)),cell_b2Lt,cell_f)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_gx[1],1,1,1) == false
@test is_zero_block(cell_gx[1],1,1,2) == true
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,1) == true
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,2) == false
r = map(i->3.0*i,cell_b2Ltx)
test_array(cell_gx,r)

cell_f = Fill(ConstantField(3),ncells)
cell_fx = lazy_map(evaluate,cell_f,cell_x)
cell_g = lazy_map(Broadcasting(Operation(*)),cell_f,cell_b2Lt)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_gx[1],1,1,1) == false
@test is_zero_block(cell_gx[1],1,1,2) == true
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,1) == true
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,2) == false
r = map(i->3.0*i,cell_b2Ltx)
test_array(cell_gx,r)

cell_g = lazy_map(Broadcasting(Operation(-)),cell_b2Lt,cell_b2Rt)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_gx[1],1,1,1) == false
@test is_zero_block(cell_gx[1],1,1,2) == false
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,1) == true
@test is_zero_block(cell_gx[1][Block(1,1,1)],1,1,2) == false
@test is_zero_block(cell_gx[1][Block(1,1,2)],1,1,1) == true
@test is_zero_block(cell_gx[1][Block(1,1,2)],1,1,2) == false
r = map((i,j)->i-j,cell_b2Ltx,cell_b2Rtx)
test_array(cell_gx,r)

cell_g = lazy_map(Broadcasting(Operation(*)),cell_b2Lt,cell_b2R)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
@test isa(cell_gx.g.value,BlockArrayCooMap)
@test is_zero_block(cell_gx[1],1,1,1) == true
@test is_zero_block(cell_gx[1],1,1,2) == true
@test is_zero_block(cell_gx[1],1,2,1) == false
@test is_zero_block(cell_gx[1],1,2,2) == true
@test is_zero_block(cell_gx[1][Block(1,2,1)],1,1,1) == true
@test is_zero_block(cell_gx[1][Block(1,2,1)],1,1,2) == true
@test is_zero_block(cell_gx[1][Block(1,2,1)],1,2,1) == true
@test is_zero_block(cell_gx[1][Block(1,2,1)],1,2,2) == false

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
