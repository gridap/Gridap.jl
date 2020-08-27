module ApplyTests

using Test
using Gridap.Arrays
using Gridap.Fields
using FillArrays
using Gridap.TensorValues

using Gridap.Fields: MockField, MockBasis, FieldOpKernel

k = FieldOpKernel(-)

v = 3.0
d = 2
f = MockField{d}(v)
g = apply_kernel_to_field(k,f)

p = 4
p1 = Point(1,2)
p2 = Point(2,1)
p3 = Point(4,3)
p4 = Point(6,1)
x = [p1,p2,p3,p4]

fx = evaluate(f,x)
∇fx = evaluate(∇(f),x)
gx = apply_kernel(k,fx)
∇gx = apply_kernel(k,∇fx)
test_field(g,x,gx,grad=∇gx)

fi = 3.0
gi = 4.5
d = 2
f = MockField{d}(fi)
g = MockField{d}(gi)
h = apply_kernel_to_field(k,f,g)
hx = apply_kernel(k,evaluate(f,x),evaluate(g,x))
∇hx = apply_kernel(k,evaluate(∇(f),x),evaluate(∇(g),x))
test_field(h,x,hx,grad=∇hx)

fi = 3.0
gi = VectorValue(4,5)
d = 2
f = MockField{d}(fi)
g = gi
h = apply_kernel_to_field(k,f,g)
hx = apply_kernel(k,evaluate(f,x),evaluate_field(g,x))
∇hx = apply_kernel(k,evaluate(∇(f),x),evaluate_field(field_gradient(g),x))
test_field(h,x,hx,grad=∇hx)

fi = 3.0
gi = [1,2,3,4,5]
d = 2
f = MockField{d}(fi)
g = gi
h = apply_kernel_to_field(k,f,g)
hx = apply_kernel(k,evaluate(f,x),evaluate_field(g,x))
∇hx = apply_kernel(k,evaluate(∇(f),x),evaluate_field(field_gradient(g),x))
test_field(h,x,hx,grad=∇hx)

end # module
