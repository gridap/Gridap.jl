module FieldGradientTests

using Gridap.Mappings
using Gridap.NewFields
using Gridap.NewFields: MockField, MockBasis, OtherMockBasis

using Test
using FillArrays
using Gridap.TensorValues


v = 3.0
d = 2
f = MockField{d}(v)

np = 4
p1 = Point(1,2)
p2 = Point(2,1)
p3 = Point(4,3)
p4 = Point(6,1)
x = [p1,p2,p3,p4]

fx = evaluate(f,x)
test_field(f,(x,),fx)

∇f = gradient(f)
∇fx = evaluate(∇f,x)
test_field(∇f,(x,),∇fx)

test_field(∇f,x,∇fx)
test_field(f,x,fx,grad=∇fx)

af = [f,f,f]
@test evaluate(af,x) == hcat(fx,fx,fx)

∇af = ∇(af)
evaluate(∇af,x)

_∇af = map(∇,af)
@test ∇(af) == _∇af
@test evaluate(_∇af,x) == evaluate(∇(af),x)


_∇∇af = map(∇,∇(af))
_∇∇afx = evaluate(_∇∇af,x)
@test ∇(∇(af)) == _∇∇af
@test evaluate(_∇∇af,x) == evaluate(∇(∇(af)),x)

f = ConstantField(v)
# f1 = ConstantField(v)
# f = [f1,f1,f1]

fx = evaluate(f,x)
test_field(f,(x,),fx)

∇f = gradient(f)
return_cache(∇f,x)
∇fx = evaluate(∇f,x)
test_field(∇f,(x,),∇fx)

∇∇f = gradient(gradient(f))

∇∇fx = evaluate(∇∇f,x)
∇∇fx = evaluate(∇∇f,x)
test_field(∇∇f,(x,),∇∇fx)

f = ConstantField(v)
# f1 = ConstantField(v)
f = [f,f,f]

fx = evaluate(f,x)
test_field(f,(x,),fx)

∇f = gradient(f)
∇fx = evaluate(∇f,x)
test_field(∇f,(x,),∇fx)
evaluate(∇f,vcat(x,x)) == vcat(∇fx,∇fx)

∇∇f = gradient(gradient(f))
∇∇fx = evaluate(∇∇f,x)
∇∇fx = evaluate(∇∇f,x)
test_field(∇∇f,(x,),∇∇fx)

end #module
