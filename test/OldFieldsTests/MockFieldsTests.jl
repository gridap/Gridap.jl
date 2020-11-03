module MockFieldsTests

using Gridap.Fields
using Gridap.Fields: MockField, MockBasis, OtherMockBasis
using Gridap.TensorValues

np = 4
p = Point(1,2)
x = fill(p,np)

v = 3.0
d = 2
f = MockField{d}(v)
fx = fill(v,np)
test_field(f,x,fx)

∇fx = fill(VectorValue(v,0.0),np)
∇f = gradient(f)

test_field(∇f,x,∇fx)
test_field(f,x,fx,grad=∇fx)

ndof = 8
b = MockBasis{d}(v,ndof)
bx = fill(v,np,ndof)
∇bx = fill(VectorValue(v,0.0),np,ndof)
test_field(b,x,bx,grad=∇bx)

b = OtherMockBasis{d}(ndof)
bx = fill(2*p,np,ndof)
∇bx = fill(TensorValue(2.0,0.0,0.0,2.0),np,ndof)
test_field(b,x,bx,grad=∇bx)

end # module
