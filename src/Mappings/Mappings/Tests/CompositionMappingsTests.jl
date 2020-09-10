module CompositionMappingsTests

using Test
# using Gridap.Arrays
using Gridap.Mappings
using FillArrays
using Gridap.TensorValues

using Gridap.NewFields
using Gridap.NewFields: MockField, MockBasis, OtherMockBasis

function foo2(t::Vararg{T}) where T
  t
end

function foo3(t...)
  t
end

function foo(t)
  t
end

foo3(1,2)
foo3(1)


k = -

v = 3.0
d = 2
f = MockField{d}(v)
g = composition(k,f)

p = 4
p1 = Point(1,2)
p2 = Point(2,1)
p3 = Point(4,3)
p4 = Point(6,1)
x = [p1,p2,p3,p4]

fx = evaluate(f,x)
∇fx = evaluate(∇(f),x)
gx = evaluate(k,fx)
∇gx = evaluate(k,∇fx)
# santiagobadia :  I don't like the different notation for
# mapping and kernel, i.e., tuple vs non-tuple, to be fixed
test_mapping(g,(x,),gx) #,grad=∇gx)

# @santiagobadia : Create the gradient method for composition
# using chain rule
∇g = composition(k,∇(f))
test_mapping(∇g,(x,),∇gx)

fi = 3.0
gi = 4.5
d = 2
f = MockField{d}(fi)
g = MockField{d}(gi)

# @santiagobadia : I like the tuple, k ∘ (f,g)
h = composition(k,f,g)

hx = evaluate(k,evaluate(f,x),evaluate(g,x))
∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
test_mapping(h,(x,),hx) #,grad=∇hx)

∇h = composition(k,∇(f),∇(g))
test_mapping(∇h,(x,),∇hx)

fi = 3.0
gi = VectorValue(4,5)
d = 2
f = MockField{d}(fi)
g = ConstantField(gi)
h = composition(k,f,g)

hx = evaluate(k,evaluate(f,x),evaluate(g,x))
∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
test_mapping(h,(x,),hx)

∇h = composition(k,∇(f),∇(g))
test_mapping(∇h,(x,),∇hx)


k = BroadcastMapping(-)
fi = 3.0
gi = [1,2,3,4,5]
d = 2
f = MockField{d}(fi)
g = ConstantField(gi)
h = composition(k,f,g)

evaluate(f,x)
evaluate(g,x)

# @santiagobadia : A field can be paralelised wrt array of points
# A mapping on the contrary can take multiple arguments...
# We must create two evalutes, one for arrays of tuples and
# one for one tuple
# @santiagobadia : The previous version of gradient is wrong
# @santiagobadia : We could also have derivative (transpose of gradient),
# i.e., a Jacobian for f : Rn -> Rm

hx = evaluate(k,evaluate(f,x),evaluate(g,x))
test_mapping(h,(x,),hx) #,grad=∇hx)

∇h = composition(k,∇(f),∇(g))
∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
test_mapping(∇h,(x,),∇hx)

end # module
