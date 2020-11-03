module ComposeTests

using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: MockField
using Gridap.Fields: Comp
using FillArrays

import Gridap.Fields: ∇

np = 4
p = Point(1,2)
x = fill(p,np)

fun(x,y) = 2*x
∇fun(x,y) = VectorValue(2*one(x[1]),2*one(x[1]))
∇(::typeof(fun)) = ∇fun

v = 3.0
d = 2
f = MockField{d}(v)
fx = fill(v,np)

k = Comp(fun)
test_mapping(k,(v,v),2*v)

g = compose(fun,f,f)
gx = 2*fx
∇gx = fill(∇fun(v,v),np)
test_field(g,x,gx,grad=∇gx)

l = 10
af = Fill(f,l)
ax = fill(x,l)
ag = compose(fun,af,af)
agx = fill(gx,l)
a∇gx = fill(∇gx,l)
test_array_of_fields(ag,ax,agx,grad=a∇gx)

np = 4
p = Point(1,2)
x = fill(p,np)

v = VectorValue(1,2,3)
d = 2
f = MockField{d}(v)
w = 4
d = 3
g = MockField{d}(w)

fx = evaluate(f,x)
gx = evaluate(g,fx)
∇gx = evaluate(∇(g),fx)

h = compose_fields(g,f)
test_field(h,x,gx,grad=∇gx)

h = compose(g,f)
test_field(h,x,gx,grad=∇gx)

l = 10
af = Fill(f,l)
ag = Fill(g,l)
ax = Fill(x,l)

afx = evaluate(af,ax)
agx = evaluate(ag,afx)
∇agx = evaluate(∇(ag),afx)

ah = compose_field_arrays(ag,af)
test_array_of_fields(ah,ax,agx,grad=∇agx)

ah = compose(ag,af)
test_array_of_fields(ah,ax,agx,grad=∇agx)

end # module
