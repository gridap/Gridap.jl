module CompositionMappingsTests

using Test
# using Gridap.Arrays
using Gridap.Mappings
using FillArrays
using Gridap.TensorValues

using Gridap.NewFields
using Gridap.NewFields: MockField, MockBasis, OtherMockBasis

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
fi = VectorValue(3.0,0.0)
d = 2
f = MockField{d}(fi)
g = broadcast(ConstantField,x)
h = composition(k,f,g)

# @santiagobadia : A field can be paralelised wrt array of points
# A mapping on the contrary can take multiple arguments...
# We must create two evaluates, one for arrays of tuples and
# one for one tuple
# @santiagobadia : The previous version of gradient is wrong
# @santiagobadia : We could also have derivative (transpose of gradient),
# i.e., a Jacobian for f : Rn -> Rm

hx = evaluate(k,evaluate(f,x),evaluate(g,x))
test_mapping(h,(x,),hx) #,grad=∇hx)

G = ∇(g)

v = 3.0
v = VectorValue(3.0,0.0)
f = MockField{d}(v)

f = ConstantField(v)

fx = evaluate(f,x)
test_field(f,(x,),fx)

# cf =
# @enter return_cache(f,x)

# c = return_gradient_cache(f,x)
# evaluate_gradient!(c,f,x)
# evaluate_gradient!(c,f,vcat(x,x))


# ∇f = gradient(f)



∇fx = evaluate(∇f,x)
test_field(∇f,(x,),∇fx)










f = ∇f
w = evaluate(f,x)
v = ∇fx
np, = size(w)
@test length(x) == np
@test ==(w,v)
@test typeof(w) == return_type(f,x)

cf = return_cache(f,x)
cf
r = evaluate!(cf,f,x)
@test ==(r,v)

_x = vcat(x,x)
_v = vcat(v,v)
_w = evaluate!(cf,f,_x)
@test ==(_w,_v)

_w .== _v
size(_w)

_v
size(_v)
∇(f)

evaluate(G,x)

evaluate(∇(G),x)

H = ∇(∇(g[1]))
return_cache(H,x)


return_cache(H,x)
return_hessian_cache(g,x)

evaluate(∇(g[1]),x)
evaluate(∇(∇(g[1])),x)












∇∇f = ∇(∇(f))
∇∇fx = evaluate(∇∇f,x)
test_field(∇∇f,x,∇∇fx)

test_field(f,x,fx,grad=∇fx,hessian=∇∇fx)






∇h = composition(k,∇(f),∇(g))
∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
test_mapping(∇h,(x,),∇hx)

end # module
