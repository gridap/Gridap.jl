module FieldArraysTests

using Test
using Gridap.Arrays
using Gridap.Mappings
using Gridap.NewFields
using FillArrays
using Gridap.TensorValues


np = 4
p = Point(1,2)
x = fill(p,np)

v = 3.0
d = 2
f = MockField{d}(v)
fx = fill(v,np)
∇fx = fill(VectorValue(v,0.0),np)
test_field(f,x,fx,grad=∇fx)

nf = 3
b = Fill(f,nf)
∇b = gradient(b)
∇bx = evaluate(∇b,x)
bx = fill(v,np,nf)
∇bx = fill(VectorValue(v,0.0),np,nf)
test_field(b,x,bx,grad=∇bx)

l = 10
af = Fill(f,l)
ax = fill(x,l)
afx = fill(fx,l)
a∇fx = fill(∇fx,l)
test_mapped_array(af,ax,afx,grad=a∇fx)

l = 10
ab = Fill(b,l)
ax = fill(x,l)
abx = fill(bx,l)
a∇bx = fill(∇bx,l)
test_mapped_array(ab,ax,abx,grad=a∇bx)

# s = Fill(BroadcastMapping(+),l)
s = Fill(+,l)
ag = apply_mapping(composition,s,af,af)
gx = fill(v+v,np)
∇gx = fill(VectorValue(v+v,0.0),np)
agx = fill(gx,l)
a∇gx = fill(∇gx,l)
# test_mapped_array(ag,ax,agx,grad=a∇gx)
# @santiagobadia : Still need to define the gradient of the
# composition mapping... it is very operator dependent, it
# must be implemented in each case. It is obvious to implement
# what we have now, but this is mathematically wrong
test_mapped_array(ag,ax,agx)

∇af = apply_mapping(gradient,af)
∇ag = apply_mapping(composition,s,∇af,∇af)
test_mapped_array(∇ag,ax,a∇gx)

# @santiagobadia : Not sure what we are testing here
# struct FieldPlaceHolder <: NewField end
# ag = apply_to_field_array(FieldPlaceHolder,bcast(+),af,af)
# test_array(evaluate_field_array(ag,ax),agx)

w = 2.0
aw = fill(ConstantField(w),l)
ag = apply_mapping(composition,s,af,aw)
gx = fill(v+w,np)
∇gx = fill(VectorValue(v,0.0),np)
agx = fill(gx,l)
a∇gx = fill(∇gx,l)
test_mapped_array(ag,ax,agx)#,grad=a∇gx)

∇af = apply_mapping(gradient,af)
∇aw = apply_mapping(gradient,aw)
∇ag = apply_mapping(composition,s,∇af,∇aw)
test_mapped_array(∇ag,ax,a∇gx)

# ag = apply_to_field_array(FieldPlaceHolder,bcast(+),af,aw)
# test_array(evaluate_field_array(ag,ax),agx)

l = 10
af = Fill(f,l)
ax = Fill(x,l)
s = Fill(+,l)
ag = apply_mapping(composition,s,af,af)
r1 = apply_mapping(ag,ax)
@test isa(r1,Fill)

np = 4
p = Point(1,2)
x = fill(p,np)
v = 2.0
d = 2
ndof = 8
wi = 3.0
w = fill(ConstantField(wi),ndof)
r = fill(v+wi,np,ndof)
f = MockBasis(d,v,ndof)
∇fx = evaluate(∇(f),x)
af = Fill(f,l)
ax = fill(x,l)
aw = fill(w,l)
s = Fill(+,l)
ag = apply_mapping(composition,s,af,aw)
agx = fill(r,l)

# @santiagobadia : Next step, arrays of fields working
# test_mapped_array(ag,ax,agx)#,grad=a∇gx)


# evaluate(ag[1],ax[1])

# typeof(ag[1])

# afx = apply_mapping(af,ax)
# awx = apply_mapping(aw,ax)

# a∇gx = fill(∇fx,l)
# ∇af = apply_mapping(gradient,af)
# ∇aw = apply_mapping(gradient,aw)
# ∇ag = apply_mapping(composition,s,∇af,∇aw)
# test_mapped_array(∇ag,ax,a∇gx)

v = 2.0
d = 2
wi = 3.0
# w = fill(wi,ndof)
w = fill(ConstantField(wi),ndof)
r = fill(v+wi,np,ndof)
f = MockField{d}(v)
∇r = fill(VectorValue(v,0.0),np,ndof)
af = Fill(f,l)
ax = fill(x,l)
aw = fill(w,l)
s = Fill(BroadcastMapping(+),l)
ag = apply_mapping(composition,s,af,aw)
agx = fill(r,l)
# @santiagobadia : To be fixed
# test_mapped_array(ag,ax,agx)#,grad=a∇gx)

typeof(af)
a∇gx = fill(∇r,l)
∇af = apply_mapping(gradient,af)
# @santiagobadia : To be fixed
# ∇aw = apply_mapping(gradient,aw)
# ∇ag = apply_mapping(composition,s,∇af,∇aw)
# test_mapped_array(∇ag,ax,a∇gx)

# lazy_append

np = 4
p = Point(1,2)
x = fill(p,np)

v = 3.0
d = 2
f = MockField{d}(v)
fx = evaluate(f,x)
∇fx = evaluate(∇(f),x)

l = 10
af = Fill(f,l)
ax = fill(x,l)
afx = fill(fx,l)
a∇fx = fill(∇fx,l)

np = 5
p = Point(4,3)
x = fill(p,np)

v = 5.0
d = 2
f = MockField{d}(v)
fx = evaluate(f,x)
∇fx = evaluate(∇(f),x)

l = 15
bf = Fill(f,l)
bx = fill(x,l)
bfx = fill(fx,l)
b∇fx = fill(∇fx,l)

cf = lazy_append(af,bf)
cx = lazy_append(ax,bx)

apply_mapping(af,ax)

cfx = apply_mapping(cf,cx)
∇cf = apply_mapping(gradient,cf)
∇cfx = apply_mapping(∇cf,cx)
rfx = vcat(afx,bfx)
∇rfx = vcat(a∇fx,b∇fx)
test_mapped_array(cf,cx,rfx)#,grad=∇rfx)
test_mapped_array(∇cf,cx,∇rfx)#,grad=∇rfx)

# @santiagobadia : Do we need this?
# @test isa(cfx, AppendedArray)
# @test isa(∇cfx, AppendedArray)

end # module
