module LinearCombinationTests

using Test
using Gridap.Arrays
using Gridap.NewKernels
using Gridap.NewFields
using Gridap.NewFields: MockField, MockBasis
using FillArrays
using Gridap.TensorValues



a = subtypes(Kernel)
a[1:10]

np = 4
p = Point(1,2)
x = fill(p,np)
v = 2.0
d = 2
ndof = 8
wi = 3.0
w = fill(wi,ndof)
f = MockBasis{d}(v,ndof)
g = linear_combination(f,w)
fx = evaluate(f,x)
∇fx = evaluate(∇(f),x)
gx = fx*w
∇gx = ∇fx*w
test_field(g,x,gx,grad=∇gx)

l = 10
af = Fill(f,l)
ax = fill(x,l)
aw = fill(w,l)
ag = lincomb(af,aw)
agx = fill(gx,l)
a∇gx = fill(∇gx,l)
test_array_of_fields(ag,ax,agx,grad=a∇gx)

l = 0
af = Fill(f,l)
ax = fill(x,l)
aw = fill(w,l)
ag = lincomb(af,aw)
agx = fill(gx,l)
a∇gx = fill(∇gx,l)
test_array_of_fields(ag,ax,agx,grad=a∇gx)

end #module
