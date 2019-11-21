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
test_kernel(k,(v,v),2*v)

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

end # module
