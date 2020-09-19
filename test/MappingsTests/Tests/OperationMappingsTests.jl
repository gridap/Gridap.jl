module OperationMappingsTests

using Test
# using Gridap.Arrays
using Gridap.Mappings
using FillArrays
using Gridap.TensorValues

using Gridap.NewFields
using Gridap.NewFields: MockField, MockBasis, OtherMockBasis

x = [1,2]

fa(x) = 2*x
fb(x) = sqrt.(x)
fab = operation(fa,fb)
test_mapping(fab,(x,),[2.0,sqrt(2)*2.0])

p = 4
p1 = [1,2]
p2 = [2,1]
p3 = [4,3]
p4 = [6,1]

x = [p1,p2,p3,p4]

aa = Fill(fa,4)
bb = Fill(fb,4)
cm = apply(operation,aa,bb)
r = apply(cm,x)
@test all([ r[i] ≈ 2*(sqrt.(x[i])) for i in 1:4])

# k = -

# v = 3.0
# d = 2
# f = MockField{d}(v)
# g = composition(k,f)

# p = 4
# p1 = Point(1,2)
# p2 = Point(2,1)
# p3 = Point(4,3)
# p4 = Point(6,1)
# x = [p1,p2,p3,p4]

# fx = evaluate(f,x)
# ∇fx = evaluate(∇(f),x)
# gx = evaluate(k,fx)
# ∇gx = evaluate(k,∇fx)
# test_mapping(g,(x,),gx) #,grad=∇gx)

# # @santiagobadia : Create the gradient method for composition
# # using chain rule
# ∇g = composition(k,∇(f))
# test_mapping(∇g,(x,),∇gx)


# g = field_composition(k,f)
# @test g isa NewField
# gx = evaluate(g,x)
# cache = return_cache(g,x)
# test_mapping(g,(x,),gx) #grad=∇gx)

# fi = 3.0
# gi = 4.5
# d = 2
# f = MockField{d}(fi)
# g = MockField{d}(gi)

# # @santiagobadia : I like the tuple, k ∘ (f,g)
# h = composition(k,f,g)

# hx = evaluate(k,evaluate(f,x),evaluate(g,x))
# ∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
# test_mapping(h,(x,),hx) #,grad=∇hx)

# ∇h = composition(k,∇(f),∇(g))
# test_mapping(∇h,(x,),∇hx)

# fi = VectorValue(3.0,0.0)
# gi = VectorValue(4,5)
# d = 2
# f = MockField{d}(fi)
# g = ConstantField(gi)
# h = composition(k,f,g)

# hx = evaluate(k,evaluate(f,x),evaluate(g,x))
# ∇hx = evaluate(k,evaluate(∇(f),x),evaluate(∇(g),x))
# test_mapping(h,(x,),hx)

# ∇h = composition(k,∇(f),∇(g))
# test_mapping(∇h,(x,),∇hx)


# k = BroadcastMapping(-)
# fi = VectorValue(3.0,0.0)
# d = 2
# f = MockField{d}(fi)
# g = broadcast(ConstantField,x)
# h = composition(k,f,g)

# # @santiagobadia : A field can be paralelised wrt array of points
# # A mapping on the contrary can take multiple arguments...
# # We must create two evaluates, one for arrays of tuples and
# # one for one tuple
# # @santiagobadia : The previous version of gradient is wrong
# # @santiagobadia : We could also have derivative (transpose of gradient),
# # i.e., a Jacobian for f : Rn -> Rm

# hx = evaluate(k,evaluate(f,x),evaluate(g,x))
# test_mapping(h,(x,),hx) #,grad=∇hx)

# G = ∇(g)

# v = 3.0
# v = VectorValue(3.0,0.0)
# f = MockField{d}(v)

# f = ConstantField(v)

# fx = evaluate(f,x)
# test_field(f,(x,),fx)

end # module
